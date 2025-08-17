export find_jwt


function find_jwt(wtd::FT, z₋ₕ::Vector{FT}) where {FT<:Real}
  nzg = length(z₋ₕ) - 1
  jwt = 1
  for k in 1:nzg
    if wtd < z₋ₕ[k]
      jwt = k # 地下水水位所在的上一层
      break
    end
  end
  return jwt
end
