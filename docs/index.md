---
title: Overview
---
# Overview

TITAN (Treatment of Infectious Transmissions through Agent-based Network) is an agent-based simulation model used to explore contact transmission in complex social networks. Starting with the initializing agent population, TITAN iterates over a series of stochastic interactions where agents can interact with one another, transmit infections through various medium, and enter and exit the care continuum. The purpose of TITAN is to evaluate the impact of prevention and treatment models on incidence and prevalence rates of the targeted disease(s) through the use of data fitting simulated trajectories and rich statistics of primary/sub-population attributable proportions.

Agent populations are defined as graphs (nodes connected by edges). Nodes in the graph are used to represent the attributes (or collection of attributes) of an agent (person), and edges define the type of relationship between agents. In practice, a graph represents a social network of connected people through various relationship types, and provides the medium for which agents can interact.

<div>
  <iframe width="100%" height="814" frameborder="0"
    src="https://observablehq.com/embed/@mcmcgrath13/visualizing-titan?cells=viewof+circle_view"></iframe>
  <figcaption>Visualizing the TITAN network during a model run.  Hover over an agent's dot to highlight their relationships to other agents. See <a href="https://observablehq.com/@mcmcgrath13/visualizing-titan">the Visualizing TITAN Observable notebook</a> to learn more.</figcaption>
</div>
