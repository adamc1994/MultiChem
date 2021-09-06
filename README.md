# MultiChem Overview
To compare the performance of multi-objective optimisation algorithms (MOOAs) for optimising chemical systems, a kinetic-based reaction simulator has been created. Four reactions with known kinetic parameters (pre-exponential factors and activation energies) were identified in the literature: (i) Van de Vase reaction; (ii) nucleophilic aromatic substitution between 2,4-difluoronitrobenzene and morpholine; (iii) isomerisation of lactose to lactulose; (iv) Paal-Knorr reaction between 2,5-hexanedione and ethanolamine. These examples provide a good representation of non-competitive (iv) and competitive reactions, including competing parallel (i, ii & iii) and consecutive pathways (ii & iii). Although reactions (iii) and (iv) contain reversible reactions, the k-1 rate constants are negligible and are therefore omitted. Six test problems have been formulated using reaction variable limits as boundaries and different process metrics as objectives. The boundaries and objectives for each problem have been selected to ensure a diverse range of Pareto fronts were generated in terms of morphology, uniformity and continuity.
|Test Problem |	Variables |	Objectives |	Description of Pareto Front|
|-------------|------------|------------|-----------------------------|
|VdV1 |	2 |	2 |	Density of solutions fall away near to the Pareto front, which is non-uniformly distributed between linear and convex regions|
|SNAr1 |	2 |	3 |	Optimal solutions follow a convoluted path through objective space with concave regions |
|SNAr2 |	4 |	3	| Convex, non-uniformly distributed Pareto front|
|Lactose1 |	2 |	2 |	Pareto front is a convex curve with many solutions |
|PK1 |	2 |	2 |	Pareto front is a convex curve with relatively few solutions |
|PK2 |	3 |	2 |	Pareto front consists of three discontinuous linear and concave regions |
## How to use
Full instructions on how to use the package can be found in the ![README](Test%20Problems%20Code/README.docx) document

