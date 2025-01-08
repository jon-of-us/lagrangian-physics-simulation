module Pendulum

import DataStructures as ds
include("./Settings.jl")
import .Settings as s
import GLMakie as mk

export 
  State, lagrangian, normalize_energy!, 
  plot!, update_plot!


#
# Pendulum system
#
mutable struct State
  # physics
  """angles"""
  x::Vector{Float64}
  """angular velocities"""
  v::Vector{Float64}
  total_energy::Float64

  # plotting obeservables
  rods::mk.Observable{Vector{mk.Point2f0}}
  balls::mk.Observable{Vector{mk.Point2f0}}
  trace::mk.Observable{ds.CircularBuffer{mk.Point2f0}}
end

#
# physics
#
function position(θ)
  pos = [zeros(2)]
  for i in 1:s.n_rods
    rod = s.rod_lengs[i] * [cos(θ[i]), sin(θ[i])]
    push!(pos, pos[end] + rod)
  end
  return vcat(pos...)
end

function velocity(θ, dθ)
  vel = [zeros(2)]
  for i in 1:s.n_rods
    dRod = s.rod_lengs[i] * [-sin(θ[i]), cos(θ[i])]
    push!(vel, vel[end] + dRod * dθ[i])
  end
  return vcat(vel...)
end

function lagrangian(θ, dθ)
  euclid_x = position(θ)
  euclid_v = velocity(θ, dθ)
  T = 0.5 * sum(euclid_v .^ 2)
  V = s.g * sum(euclid_x[2:2:end])
  L = T - V
  return L
end

function energy(θ, dθ)
  x = position(θ)
  v = velocity(θ, dθ)
  T = 0.5 * sum(v .^ 2)
  V = s.g * sum(x[2:2:end])
  return T + V
end

function normalize_energy!(pendulum::State)
  euclid_x = position(pendulum.x)
  euclid_v = velocity(pendulum.x, pendulum.v)
  T = 0.5 * sum(euclid_v .^ 2)
  V = s.g * sum(euclid_x[2:2:end])
  should_be_T = pendulum.total_energy - V
  if T > 1e-3
    pendulum.v *= sqrt(abs(should_be_T / T))
  end
end

#
# rendering
#
function to_point(pos)
  return [
    mk.Point2f0(pos[2*i-1:2*i])
    for i in 1:s.n_rods+1
  ]
end

function new()
  x = s.start_θ
  v = zeros(s.n_rods)
  start_energy = energy(x, v)
  pos = to_point(position(x))
  trace = ds.CircularBuffer{mk.Point2f0}(s.trajectory_length)
  fill!(trace, pos[end])
  return State(
    x,
    v,
    start_energy,
    mk.Observable(pos),
    mk.Observable(pos[2:end]),
    mk.Observable(trace)
  )
end

function plot!(pendulum::State)
  fig = mk.Figure()
  display(fig)
  ax = mk.Axis(
    fig[1, 1],
    title="Pendulum",
    aspect=mk.DataAspect(),
  )
  mk.xlims!(ax, -sum(s.rod_lengs), sum(s.rod_lengs))
  mk.ylims!(ax, -sum(s.rod_lengs), sum(s.rod_lengs))
  mk.lines!(ax, pendulum.rods, color=s.line_color)
  mk.scatter!(ax, pendulum.balls, color=s.ball_color)
  mk.lines!(ax, pendulum.trace, color=s.trace_color)
end

function update_plot!(pendulum::State)
  pos = to_point(position(pendulum.x))
  pendulum.rods[] = pos
  pendulum.balls[] = pos[2:end]
  push!(pendulum.trace[], pos[end])
  pendulum.trace[] = pendulum.trace[]
end

end