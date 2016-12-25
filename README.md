The travelling salesman problem ([TSP](https://en.wikipedia.org/wiki/Travelling_salesman_problem)) GUI ([Xlib](https://en.wikipedia.org/wiki/Xlib)) simulation using the [Held–Karp algorithm](https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm) (dynamic programming)

Has two input modes:
* Command-Line: Pass a .txt file containing integers 0 to 500
with two coordinates "x, y" per line.

![alt text](https://github.com/NikaTsanka/tsp-held-karp/blob/master/demo-cmd.gif "Mode: Command-Line")

* Interactive: The coordinates can be entered by mouseclicks (left mouseclick enters a point, a right mouseclick ends the input phase).

![alt text](https://github.com/NikaTsanka/tsp-held-karp/blob/master/demo-interactive.gif "Mode: Interactive")

Use the right mouseclick to display the solution.

If the number of coordinates is at most fifteen, optimum TSP tour is computed.

If the number is greater than fifteen, the coordinates are sorted according to x axis, divided into groups of at most fifteen,
 and solved each group individually, after that they are connected left-to-right sequence optimally
 
Cool right? :octocat: 


