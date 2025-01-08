
import LinearAlgebra as la
Base.show(io::IO, f::Float64) = Printf.@printf(io, "%s", abs(f) < 1e-12 ? "0" : (1e-4 < abs(f) < 1e3 ? Printf.@sprintf("%.3f", f)[1:5] : Printf.@sprintf("%.4e", f)))

include("./Settings.jl");
import .Settings as s
include("./Pendulum.jl");
import .Pendulum as sys

function grad(x, f)
  n = length(x)
  ∇f = zeros(n)
  fx = f(x)
  for i in 1:n
    x2 = copy(x)
    x2[i] += s.δ
    ∇f[i] = (f(x2) - fx) / s.δ
  end
  return ∇f
end

function physics_step!(system)
  for _ in 1:s.n_simsteps
    x = system.x
    v = system.v
    n = length(x)

    ∇vL = grad(v, v -> sys.lagrangian(x, v))
    ∇xL = grad(x, x -> sys.lagrangian(x, v))

    step_x = x + s.δ * v 
    step_x_∇vL = grad(v, v -> sys.lagrangian(step_x, v))
    dtx_∇vL = (step_x_∇vL - ∇vL) / s.δ

    Jy_∇vL = zeros(n, n)
    for i in 1:n
      v2 = copy(v)
      v2[i] += s.δ
      ∇vL2 = grad(v2, v -> sys.lagrangian(x, v))
      Jy_∇vL[:, i] = (∇vL2 - ∇vL) / s.δ
    end 
    Jy_∇vL += rand(n, n) * s.δ * 1e-10
    dv = Jy_∇vL \ (∇xL - dtx_∇vL)

    system.x += s.simstep * v 
    system.v += s.simstep * dv
    sys.normalize_energy!(system)
  end
end

function animation_step!(system)
  physics_step!(system)
  sys.update_plot!(system)
end;

system = sys.new()
fig = sys.plot!(system)

while true
  animation_step!(system)
  sleep(s.timestep)
end
