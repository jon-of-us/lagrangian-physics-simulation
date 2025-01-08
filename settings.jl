module Settings

import GLMakie as mk

# constants
# const rod_lengs = [3, 2, 4, 3]
# const start_θ = [0, 1, 2, 1]
n = 3
const rod_lengs = rand(n) .+ 0.5
const start_θ = rand(n) 
const n_rods = length(rod_lengs)
const g = 9.81

# style 
const line_width = 2
const hinge_size = 3
const trajectory_length = 100
color(α) = mk.RGBA(0.1, 0.1, 0.1, α)
line_color = color(1)
ball_color = color(1)
trace_color = [
  color((i / trajectory_length)^2)
  for i in 1:trajectory_length
]

#computation
const fps = 30
const timestep = 1 / fps
const n_simsteps = 1
const simstep = timestep / n_simsteps
const δ = 1e-3

end
