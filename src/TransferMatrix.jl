function TransferMatrix(lens::Lens)
    (; M) = lens
    τ, ϕ = eachcol(M)
    M = prod([1.0 τ[i]; -ϕ[i] 1.0-τ[i]*ϕ[i]] for i = reverse(axes(M, 1)))
    TransferMatrix(M)
end

extend(M, τ, τ′) = [1.0 τ′; 0.0 1.0] * M * [1.0 τ; 0.0 1.0]

transfer(M::Matrix, v::Vector, τ, τ′) = extend(M, τ, τ′) * v
transfer(system::System, v::Vector, τ, τ′) = transfer(system.M, v, τ, τ′)
transfer(M::TransferMatrix, v::Vector, τ, τ′) = transfer(M.M, v, τ, τ′)

reverse_transfer(M::Matrix, v::Vector, τ′, τ) = extend(M, τ, τ′) \ v

function reverse_transfer(system::System, v::Vector, τ′, τ)
    reverse_transfer(system.M, v, τ′, τ)
end

function reverse_transfer(M::TransferMatrix, v::Vector, τ′, τ)
    reverse_transfer(M.M, v, τ′, τ)
end

flatten(transfer_matrix::TransferMatrix) = flatten(transfer_matrix.M)

function flatten(M::Matrix)
    f = -inv(M[2,1])
    EFFD = -M[2,2] * f
    EBFD = M[1,1] * f
    return (; f, EFFD, EBFD)
end
