# ChemotaxisWorld
Developing a chemotaxis worldfile for MABE project

	The simulated environment is an xy plane in which an attractant gradient increases when moving in the positive x direction. "Cells" spawn at the origin with a random orientation between 0 and 2pi radians, which is measured counter-clockwise from the positive x axis. These cells are simulated for eval_ticks iterations before having their data reported to MABE. During each tick the cell can either tumble or run. A run moves the cell forward by speed units in the direction the cell is pointing, then slightly modifies the direction by adding a rotational diffusion coefficient. A tumble sets the cell's direction randomly between [0,2pi] radians. The likelyhood of the cell entering the tumble state is the tumble bias, which is calculated from the cell’s output.

	The attractant gradient can be set to increase either linearly or exponentially. If a linear gradient is selected, then the “concentration” of attractant is calculated as conc = (x_position * slope) + base. Alternatively the exponential gradient is calculated by conc = base + e^(slope * x_position). The cell cannot measure the gradient directly. Instead, the cell has a 16-bit input that conveys the magnitude of the relative change in concentration since the previous tick. The number of ones in the 16 bit input is calculated by finding the relative change in concentration between ticks, then multiplying it by a “multiplier”. The cell can specify this multiplier, really an exponent, to be between [0,31]. The number of ones to add is ((conc – oldconc)/oldconc)*(2^multiplier). The resulting number is rounded and then added to 8. This is done to have the reading be 8 ones whenever the concentration stays the same. The number of ones decreases if the concentration drops, and the number of ones increases if the concentration goes up. The magnitude of the change, and effectively the sensitivity, is modified by the cell’s control of the multiplier. This resulting number of ones could potentially be hugely positive or negative, 2^31 is approximately 2 billion, while the number of inputs is only sixteen. So we limit the range to [0,16]. The ones are then fed into the cell’s Markov Brain, and if necessary padded with 0s to reach 16 bits. The brain is then allowed to update for a user-specified number of iterations before the output is read.

	The cells output consists of two items: the tumble bias and the multiplier. The tumble bias has 16 nodes of input and the multiplier takes 8 nodes of input. These outputs are not necessarily bits; the Markov Brain implementation in MABE does not appear to OR the inputs to a single node, instead it adds them. Both items are simply the sums of their respective nodes. This value is sanity checked to prevent huge and negative numbers from entering the equation. The multiplier is limited to integers in the range [0,31], while the tumble bias is allowed to range from [1/17,1]. The bias is simply the chance that the cell has of entering the tumbling state. 

	Once the tumble bias is known, the cell can then enter a tumble or run. Movement is calculated as expected for a 2D plane. The x displacement is dx = cos(theta)*speed, and the y displacement is dy = sin(theta)*speed. Theta is the cell’s angle from the positive x axis in radians. A tumble randomly reorients the cell. 

	The cell is scored on how far it moves in the positive x direction. Each cell should always be tested more than once to reduce noise. Runs can be noisy because the cell can spawn oriented in any direction, and there is considerable chance involved in the outcome of tumbles and the effect of rotational diffusion. This can be set in the settings_world file value “repeats”. Around ten attempts appears to be reasonable when variable_environment is not set. Even then, the most performant organism in a run will often be much better than the average performance. It is then difficult to determine how much chance was involved in their final score. 


	The world supports environmental variability, which causes the world’s parameters to fluctuate by a user-specified amount. The tunable parameters are slope, base concentration, and rotational diffusion coefficient. Speed is technically listed, but it is disabled in source code. Small floating point errors cause it to behave strangely and it can make comparisons difficult. For example, they can appear to move ‘farther than possible’ even with variability_speed set to zero. 

	Visualization mode runs a genome and writes a number of statistics as text files in the same directory as the executable. These include information about the organism at each tick, such as the position and direction, concentration, multiplier, number of ones in the input, etc. 









TODO: See if we can get output into a format 

look like detection.mat 










