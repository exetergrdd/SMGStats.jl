
#### special type for crosscorrelation between bit vectors

struct BitPacked{W<:Unsigned} <: AbstractVector{Bool}
    data::Vector{W}
    nbits::Int
end

function BitPacked{W}(nbits::Int) where {W<:Unsigned}
    nwords = cld(nbits, sizeof(W)*8)
    BitPacked{W}(zeros(W, nwords), nbits)
end

function BitPacked{UInt64}(bt::BitArray{1})
    BitPacked{UInt64}(bt.chunks, length(bt))
end

function BitPacked{W}(bt::BitArray{1}) where {W}
    nwords = cld(length(bt), sizeof(W)*8)
    bp = BitPacked{W}(Vector{W}(undef, nwords), length(bt))
    
    for i = 1:length(bt)
        bp[i] = bt[i]
    end
    bp
end

@inline Base.size(bp::BitPacked{W}) where {W} = (bp.nbits,)
@inline Base.axes(bp::BitPacked) = (Base.OneTo(bp.nbits),)
@inline Base.eachindex(bp::BitPacked) = Base.OneTo(bp.nbits)
@inline Base.IndexStyle(::Type{<:BitPacked}) = IndexLinear()

@inline function checkbounds(bp::BitPacked, i)
    if 1 <= i <= bp.nbits
        return true
    else
        throw(BoundsError(bp, i))
    end
end


@inline function Base.getindex(bp::BitPacked{W}, i::Int) where {W}
    @boundscheck checkbounds(bp, i)
    wordbits = sizeof(W)*8
    word, bit = divrem(i - 1, wordbits)
    @inbounds v = (bp.data[word + 1] >> bit) & one(W) == one(W)
    v
end

@inline function Base.setindex!(bp::BitPacked{W}, v::Bool, i::Int) where {W}
    @boundscheck checkbounds(bp, i)

    wordbits = sizeof(W)*8
    word, bit = divrem(i - 1, wordbits)
    mask = one(W) << bit
    @inbounds if v
        bp.data[word + 1] |= mask
    else
        bp.data[word + 1] &= ~mask
    end
    return v
end

@inline totalbits(bp::BitPacked{W}) where {W} = sum(count_ones, bp.data)


@inline function StatsBase.crosscor(x::BitPacked{W}, y::BitPacked{W}, lags::AbstractVector{Int}) where {W}
    result   = zeros(Float64, length(lags))
    crosscor!(result, x, y, lags)
end

@inline function StatsBase.crosscor!(result::AbstractVector{<:Real}, x::BitPacked{W}, y::BitPacked{W}, lags::AbstractVector{<:Integer}) where {W}
    wordbits = sizeof(W)*8
    nwords   = length(x.data)
    # result   = zeros(Int, length(lags))

    @inbounds for (j, lag) in enumerate(lags)
        wordshift, bitshift = divrem(lag, wordbits)
        s = 0
        for i in 1:(nwords - wordshift)
            a = x.data[i]
            b = y.data[i + wordshift]
            if bitshift > 0 && i + wordshift < nwords
                ### build b as a combination of current word shifted to right and next word shifted in as in necessary
                b = (b >> bitshift) | (y.data[i + wordshift + 1] << (wordbits - bitshift))
            else
                b = b >> bitshift
            end
            s += count_ones(a & b)
        end
        result[j] = s
    end
    result ./= sqrt(totalbits(x)*totalbits(y))
    result
end