# Experimental Notes

## Goel

### Notes on implementation

#### removeLoop
- Note that the implementation might affect the overall complexity.

### Tests
1. Randomized input of different sizes
2. Adversarial input--graphs that already have a large matching, but far from perfect

## Dani-Hayes

### Notes on implementation

### Tests
1. Usual scaling tests
2. Non-regular graph: a star and a chain.
3. Thm 6.1: two cliques.
4. A graph that is already 3/4 matched.

**Observation:** In case 3 and 4, 
- If we let the algorithm solve for an approximate matching, 
then it worked well. 
- If we force it to solve for a perfect matching, the algorithm still 
finished the task of test 3 comparably to test 1. But in test 4, 
completeing the 75% matching took significantly longer than starting from stratch.