# ARCM: Multi-Objective Optimization for Integrated Energy Systems

Description: A day-ahead dispatch optimization tool for Integrated Energy Systems. This project implements the Alternating Refined Constraint Method (ARCM) to provide robust, quantitative trade-off analysis between operational expenditures and carbon footprints.

Key features include:

*A two-stage solution mechanism that mathematically guarantees strong Pareto optimality for all generated solutions.
*An alternating constraint strategy that applies budget constraints to each objective in turn, ensuring a more uniform distribution of Pareto points.
*High engineering interpretability that aligns with budget-driven decision-making and supports detailed "what-if" scenario analysis.

This code has strong engineering interpretability. For further explanation, refer to this paper.

The optimisation problem is solved via the toolbox YALMIP, you can find instructions for how to install it here. You will also need to install some external MILP/MIQP solver like Gurobi or Mosek, both of which have academic licenses available. Remember to also install the Matlab functionalities of that solver.

YALMIP is very easy to use. Once installed, you can learn how to use it with file "simple_EconomicDispatch_YALMIP.m", contained in this repository, and through the documentation available in YALMIP's website.

This code has been tested with MATLAB version 2023b.

If you use this code for your own work, please cite this paper:

<ol>
<li>  Yingzhao Tang, Shixing Ding, Zhigang Lu, Wei Gu, Chongqing Kang, Yijun Xu, and Ershun Du, "<b>Alternating Refined Constraint Method: A Multi-Objective Solution Approach for Energy System Optimization</b>," <i>IEEE Transactions on Power Systems</i>, <i>Eaely Eexess</i>, 2026, DOI: 10.1109/TPWRS.2026.3682532.
</ol>
