# MEAM620-Project1
Controller of Quadrotor. Path Planning. Trajectory tracking

1. The 'trajectory generator' function is used to smooth the path generated by path planning and do minimum snap interpolation.

2. The smooth part uses the function 'detectcollision', which is a function used in MEAM 520.

3. The interpolation let the drone move without stopping at the end of each segment.