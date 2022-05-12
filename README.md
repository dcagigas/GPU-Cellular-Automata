# GPU-Cellular-Automata

This code is related to paper:

Daniel Cagigas-Muñiz, Fernando Diaz-del-Rio, Jose Luis Sevillano-Ramos, Jose-Luis Guisado-Lizar, "Efficient simulation execution of cellular automata on GPU", Simulation Modelling Practice and Theory, Volume 118, 2022, ISSN 1569-190X.

This repository contains Cuda source code of several (bi-dimensional) Cellular Automata:

- Game of Life (GoL): original cellular automaton proposed by Conway.
- Forest Fire: simulation of how a fire is propagated
- Cyclic Cellular Automaton: based on David Griffeath work. In this version there are 15 states and Von Neumann distance is considered.
- WireWorld: cellular automaton to simulate electronic circuits.

Every cellular automaton has the Cuda baseline version and other special Cuda versions that offer better performance: 
- "Look-Up Table" versions that use an array to code rules: GoL, Forest Fire, Cyclic  Cellular Automaton, WireWorld.
- AN5D versions that implements the temporal blocking technique: GoL, Forest Fire and WireWorld. 
- "Packet coding" technique that codes several cells in a 'supercell': GoL, Cyclic Cellular Automaton, WireWorld.

The WireWorld cellular automaton has also a compressed version that codes two 4 bit cells into one 8 bits cell. The cellular automaton needs half of the baseline version space in memory. 
