# Agents

The goal is to implement a spatial version of the ATN (allometric trophic network) model on a grid (cellular automaton model). To do this, we take inspiration from the Ryser paper in the context folder.

For now, two things must be different from the Ryser paper. First, dispersal will initially be density-independant. Second, dispersal will only occur between the four adjacent cells. 

## Steps
1 - Understand the Ryser paper
2 - Code simple functions to implement basic diffusion on a grid. The function will take as an input a matrix of biomass, and will output the net immigration - emigration
3 - Execute a simple simulation and output snapshots of the dispersal process for one species and validate that dispersal occurs as would be expected
4 - Improve the dispersal function to also include dispersal triggers, as is done in the Ryser paper
5 - Execute a simple simulation and output snapshots of the dispersal process for one species and validate that dispersal occurs as would be expected
6 - Compare the simple dispersal with the improved dispersal from the Ryser paper

## File architecture
|-context : papers necessary to implement the spatial ATN model
|-functions : function to implement the simulation
|-thinking : documentation
|- scripts : scripts to execute simulation
|-data : external data to calibrate the model (not needed for now)
|-figures : visual outputs of the simulation model


## Constraints
Re-read AGENTS.md and tell me what is still ambiguous and what you would need to know before answering and start coding.

All code must be written in R

Always answer in English

Always confirm with me before adding a file, modifying a file or executing code

Output your thinking process into a specific file that you update in the thinking folder. Be synthetic and make sure it is easily readable and understable for graduate ecologists.
