# riddler_dice
riddler puzzle on rolling dice


The riddler puzzle of Mar 27 2020 asks the following:
Given a 6 sided die, roll it 6 times and record the number on a each side of a new blank die. Take the resulting die and repeat until all 6 sides show the same number. What is the expeced number of rolls? Generalize to N sided die.

We can view this as a probabalistic state machine, where each die configuration is a state, and rolling it N times can give a new configuration. However, there are power(N,N) configurations, so this is cumbersome  since even with N=6 there are 46656 configurations and the state machine would have 2,176,782,336‬ possible transitions.


To reduce this complexity consider states that consist of multiple configurations that are equivalent for the purposes of this puzzle.
A first way to do this is to consider states where each state comprises a number of configurations with the same number of sides that have the same value. For example, 123456 and 123465 each have 1 of each value. We can encode such states as n1 n2 ... nN representing the number of each value. Therefore we start in state 111111 and want to end up in any of the states 600000, 060000, ... , 000006. This is more efficient but there are still 462 states for the 6 sided puzzle and it grows rapidly.


A more powerful representation is simply the number of identical values, expressed in decreasing order, without any regard to the values. That is, the representations n0 n1 ... nN means there are n0 of some value and n1 of another different value, etc. We requires ni >= n[i+1] to reduce the number of states. For 6 sided there are only 11 states 111111, 21111, 2211, 222, 3111, 321, 33, 411, 42, 51, 6. The 111111 state represents 1 of each value, 321 reprsents 3 of some value, 2 of another, and 1 of the last (eg. any permuation of aaabbc for any 3 distinct values of a, b, c). This represetation is also known as the partitions of N.
Now we can associate a value d[i] with any state i representing the average number of throws to reach state i. d[last_state] = 0 and the distance for the other states can be expressed as d[i] = sum (p [i, j] * d[j]) + 1, where p [i, j] is the probabiltiy of making a transition to state j given that the die is in state i. This represents a matrix of equations x = Ax + b which we simply rearrange into (A-I) * x = -b and solve.
The next complexity is computing p[i, j]. Consider a simple example, p [411, 321]. We can move into 321 if we roll any of 3 a, 2, b, 1 c, or 3 a, 1 b, 2c, etc. In other words, any permuation of 321 rolls of a die that has 411 values. The probabilty of rolling a is 4/6, b is 1/6, and c is 1/6. So the probability of 3a,2b,1c is comb(6,4,1,1) * (4/6)^3 * (1/6)^1 * (1/6)^1. Enumerate all of the permutations and sum to get the answer. There is some trickery in getting the permutations right, for example going to state 222 has only 1 unique case, 2a, 2b, and 2c.
The code does this both in floating point for simplicity and arbitrary precision rational numbers to see if there is any observable pattern. I can't find one.

See pall.out for the answers.
