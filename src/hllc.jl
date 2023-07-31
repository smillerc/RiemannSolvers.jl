
struct HLL <: AbstractRiemannSolver end
struct HLLC <: AbstractRiemannSolver end
struct HLLCM <: AbstractRiemannSolver end


@inline function rho_ave(L::NTuple{2,T}, R) where {T}
  L1, L2 = L
  R1, R2 = R

  return (L1 * √L2 + R1 * √R2) / (√L2 + √R2)
end

function Ustar(ρ::T, u, v, w, p, E, S, Sstar) where {T}
  # Eq. 10.73 in Toro
  return ((S - u) / (S - Sstar)) * SVector{4,T}(
    ρ,
    ρ * Sstar,
    ρ * v,
    ρ * w,
    E + (Sstar - u) * (ρ * Sstar + p / (S - u))
  )
end


function flux(RS::HLLC, L::SVector{4,T}, R::SVector{4,T}, n̂::SVector{2,T}) where {T}
  ρL, uL, vL, pL = L
  ρR, uR, vR, pR = R


  sqrt_ρL = sqrt(ρL)
  sqrt_ρR = sqrt(ρR)

  û = (uL * sqrt_ρL + uR * sqrt_ρR) / (sqrt_ρL + sqrt_ρR)

  ĉ² = (
    ((cL^2 * sqrt_ρL + cR^2 * sqrt_ρR) / (sqrt_ρL + sqrt_ρR))
    +
    0.5 * (uR - uL)^2 * (sqrt_ρL * sqrt_ρR) / (sqrt_ρL + sqrt_ρR)^2
  )

  Sstar = (
    (pR - pPL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)) /
    (ρL * (SL - uL) - ρR * (SR - uR))
  )

  SL = min(uL - cL, û - ĉ)
  SR = max(uR + cR, û + ĉ)


  if 0 <= SL
  elseif SL <= 0 <= Sstar
  elseif Sstat <= 0 <= SR
    return @SVector [FL]
  elseif 0 >= SR
  end

  return nothing
end

@inline Dstar(::SVector{3}, Sstar::T) where {T} = SVector{3,T}(0, 1, Sstar)
@inline Dstar(::SVector{4}, Sstar::T) where {T} = SVector{4,T}(0, 1, 0, Sstar)
@inline Dstar(::SVector{5}, Sstar::T) where {T} = SVector{5,T}(0, 1, 0, 0, Sstar)

@inline function F(ρ, u, p, E)
  return @SVector [
    ρ * u,
    ρ * u^2 + p,
    u * (E + p),
  ]
end

@inline function F(ρ, u, v, p, E)
  return @SVector [
    ρ * u,
    ρ * u^2 + p,
    ρ * u * v,
    u * (E + p),
  ]
end

@inline function F(ρ, u, v, w, p, E)
  return @SVector [
    ρ * u,
    ρ * u^2 + p,
    ρ * u * v,
    ρ * u * w,
    u * (E + p),
  ]
end

function HLLCFlux(Sstar, S, U⃗, F⃗)
  # Eq. 10.76 in Toro
  PLR = 0.5 * (pL + pR + ρL * (SL − uL) * (Sstar − uL) + ρR * (SR − uR) * (Sstar − uR))

  # Eq. 10.75 in Toro
  return Sstar .* (S .* U⃗ - F⃗) + S * PLR * Dstar(U⃗, Sstar)
end
