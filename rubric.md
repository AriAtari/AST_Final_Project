# Grading rubric: 200 pts

## Code (45 pts)

### Code Style (15 pts)

Each item is graded on a 3-point scale:  
**3** Standard is consistently followed  
**2** Standard is mostly followed, but there are gaps  
**1** Standard is not followed; no consistency

1. Source code lines are < 80 characters in length; long-lines are broken
2. Blank lines are used to enhance readability (but are not overdone)
3. Code is structured clearly; deeply nested loops and lengthy functions are avoided, where possible
4. Variable and function names are clear, and a consistent style, e.g., "lower\_with\_under" is followed. Constants have a distinguishing style.
5. All functions have a descriptive docstring: arguments and return values for each function are clearly specified.

### Modularity (15 pts)

Each item is worth 5 points:  
**5** Consistent application of standard  
**4** Intermediate between criteria above and below  
**3** Standard is mostly followed, but not consistently  
**2** Intermediate between criteria above and below  
**1** Standard is not followed  

1. Common tasks are put into reusable functions, rather than duplicating code.
2. Functions in a module logically belong together.
3. Reaction rate can be easily modified without rewriting code.

### Efficiency (15 pts)

Each item is worth 5 points, following the scale in the section "Modularity"

1. `numpy` arrays and functions are used for numerical calculation, as opposed to inefficiently looping over python `lists`.
2. Iterators are preferred to `while` loops.
3. Numbers are handled with scientific notation, rather than as m*10**e.

## Calculation (80 pts)

### Methods (40 pts)

**5 pts** for each item

1. Temperature and density
2. Mean molecular weight
3. Heating rate
4. Differential equations
5. Central values
6. Integration loop
7. Computation of surface luminosity, including interpolation of *T*<sub>eff</sub>
8. Determination of radius

### Testing (40 pts)

**5 pts** for each item

1. Unit tests implemented (**3 pts**) and passed (**2 pts**)
    1. const
    2. eos
    3. reactions
    4. structure
    5. zams
2. Convergence testing
    1. &delta; is appropriate&mdash;converged, but not unnecessarily small
    2. &eta; is appropriate
    3. &xi; is appropriate

## Results (30 pts)

**5 pts** for each item

1. Verification that central pressure, density is 0.77 GM<sup>2</sup>/R<sup>4</sup> and 5.99 3M/(4&pi;R<sup>3</sup>)
2. Plot of HR diagram
3. Plot of central temperature against central density
4. Plot of L(r): should find that L(r) = 0.9L(R) around r/R = 0.4
5. Plot of L(m): should find that L(m) = 0.9L(M) around m/M = 0.28
6. With enhanced pp-rate: radius is about 3 times larger, L about 10 times larger.

## Report (15 pts)

1. Prose is clear and at an appropriate level (**3 pts**)
2. Report discusses all unit tests was tested (**3 pts**)
3. Discussion of convergence: Rationale for values of &delta;, &eta;, and &xi; are provided (**3 pts**)
4. Comparison with Paxton et al. (**3 pts**)
5. Discussion of consequences if &epsilon;<sub>pp</sub> were enhanced. (**3 pts**)

## Collaborative Success (30 pts)

1. Based on peer assessment survey (to be released later)
