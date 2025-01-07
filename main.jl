
# import Pkg; Pkg.add("DataStructures")
import GLMakie as mk
import DataStructures as ds
import LinearAlgebra as la
Base.show(io::IO, f::Float64) = Printf.@printf(io, "%s", abs(f) < 1e-12 ? "0" : (1e-4 < abs(f) < 1e3 ? Printf.@sprintf("%.3f", f)[1:5] : Printf.@sprintf("%.4e", f)))

include("./settings.jl");

# %%

mutable struct Pendulum
  θ::Vector{Float64}
  dθ::Vector{Float64}
  total_energy::Float64
end

function position(θ)
  pos = [zeros(2)]
  for i in 1:n_rods
    rod = rod_lengs[i] * [cos(θ[i]), sin(θ[i])]
    push!(pos, pos[end] + rod)
  end
  return vcat(pos...)
end

function velocity(θ, dθ)
  vel = [zeros(2)]
  for i in 1:n_rods
    dRod = rod_lengs[i] * [-sin(θ[i]), cos(θ[i])]
    push!(vel, vel[end] + dRod * dθ[i])
  end
  return vcat(vel...)
end

struct Observables
  rods::mk.Observable{Vector{mk.Point2f0}}
  balls::mk.Observable{Vector{mk.Point2f0}}
  trace::mk.Observable{ds.CircularBuffer{mk.Point2f0}}
end


function to_point(pos)
  return [
    mk.Point2f0(pos[2*i-1:2*i])
    for i in 1:n_rods+1
  ]
end

function lagrangian(θ, dθ)
  x = position(θ)
  v = velocity(θ, dθ)
  T = 0.5 * sum(v.^2)
  V = g * sum(x[2:2:end])
  L = T - V
  return L
end

function energy(pendulum::Pendulum)
  x = position(pendulum.θ)
  v = velocity(pendulum.θ, pendulum.dθ)
  T = 0.5 * sum(v.^2)
  V = g * sum(x[2:2:end])
  return T + V
end

function normalize_energy!(pendulum::Pendulum)
  x = position(pendulum.θ)
  v = velocity(pendulum.θ, pendulum.dθ)
  T = 0.5 * sum(v.^2)
  V = g * sum(x[2:2:end])
  should_be_T = pendulum.total_energy - V
  if T > 1e-3
    pendulum.dθ *= sqrt(abs(should_be_T / T))
  end
end

function grad(x, f)
  n = length(x)
  ∇f = zeros(n)
  fx = f(x)
  for i in 1:n
    x2 = copy(x)
    x2[i] += δ
    ∇f[i] = (f(x2) - fx) / δ
  end
  return ∇f
end

function physics_step!(pendulum::Pendulum)
  for _ in 1:n_simsteps
    x = pendulum.θ
    v = pendulum.dθ
    n = length(x)

    ∇vL = grad(v, v -> lagrangian(x, v))
    ∇xL = grad(x, x -> lagrangian(x, v))

    step_x = x + δ * v 
    step_x_∇vL = grad(v, v -> lagrangian(step_x, v))
    dtx_∇vL = (step_x_∇vL - ∇vL) / δ

    Dy_∇vL = zeros(n, n)
    for i in 1:n
      v2 = copy(v)
      v2[i] += δ
      ∇vL2 = grad(v2, v -> lagrangian(x, v))
      Dy_∇vL[:, i] = (∇vL2 - ∇vL) / δ
    end 
    Dy_∇vL += rand(n, n) * δ * 1e-10
    dv = Dy_∇vL \ (∇xL - dtx_∇vL)

    pendulum.θ += simstep * v 
    pendulum.dθ += simstep * dv
    normalize_energy!(pendulum)
  end
end

function animation_step!(pendulum::Pendulum, obs::Observables)
  physics_step!(pendulum)
  pos = to_point(position(pendulum.θ))
  obs.rods[] = pos
  obs.balls[] = pos[2:end]
  push!(obs.trace[], pos[end])
  obs.trace[] = obs.trace[]
end;


# %%
# init state
#
pendulum = Pendulum(start_θ, zeros(n_rods), 0)
pendulum.total_energy = energy(pendulum)
pos = to_point(position(pendulum.θ))
trace = ds.CircularBuffer{mk.Point2f0}(trajectory_length)
fill!(trace, pos[end])
obs = Observables(
  mk.Observable(pos),
  mk.Observable(pos[2:end]),
  mk.Observable(trace)
);

# %%
# plotting
#
fig = mk.Figure();
display(fig);
ax = mk.Axis(
  fig[1, 1],
  title="Pendulum",
  aspect=mk.DataAspect(),
)
mk.xlims!(ax, -sum(rod_lengs), sum(rod_lengs))
mk.ylims!(ax, -sum(rod_lengs), sum(rod_lengs))
mk.lines!(ax, obs.rods, color=color(1))
mk.scatter!(ax, obs.balls, color=ball_color)
mk.lines!(ax, obs.trace, color=trace_color)


while true
  animation_step!(pendulum, obs)
  sleep(timestep)
  break
end
