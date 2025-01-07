import GLMakie as mk

# constants
const rod_lengs = [3, 1, 4]
const start_θ = [0, 1, 2]
const n_rods = length(rod_lengs)
const g = 9.81

#computation
const fps = 20
const timestep = 1 / fps
const n_simsteps = 20
const simstep = timestep / n_simsteps
const δ = 1e-3

# style 
const line_width = 2
const hinge_size = 3
const trajectory_length = 100
color(α) = mk.RGBA(0.1, 0.1, 0.1, α)
ball_color = color(1)
trace_color = [
  color((i / trajectory_length)^2)
  for i in 1:trajectory_length
]

;