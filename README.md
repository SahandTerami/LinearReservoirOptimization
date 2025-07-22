# LinearReservoirOptimization

The goal of this paper is to develop an optimization framework for tuning the topology of a linear reservoir computer, making it competitive with its nonlinear counterpart while enhancing the interpretability of reservoir computing systems.

The governing equation for a conventional decoupled reservoir computer is: (See the paper for details on how it is computed based on a conventional coupled reservoir computer)
<img width="514" height="59" alt="image" src="https://github.com/user-attachments/assets/14e969a4-4d47-43b1-af0b-dcba576b3c95" />
where q is the decoupled reservoir state, Λ is a diagonal matrix of eigenvalues, c is the input weight, and γ is a positive constant.

With the input and output frequencies of: 
<img width="338" height="171" alt="image" src="https://github.com/user-attachments/assets/02addc7a-e2f1-4515-8767-d5b82a00df89" />
where a and b are the amplitudes of input and output signals, respectively. ω is the frequency of the input signal. ϕ is the phase shift in the output signal.

By solving Eq. 9,  the reservoir's state q can be expressed as in Eq. 14.
<img width="555" height="84" alt="image" src="https://github.com/user-attachments/assets/e618e284-1cfa-46c8-8257-78d4156321df" />
M and θ are according to Eq. 15 
<img width="568" height="144" alt="image" src="https://github.com/user-attachments/assets/b756e6f1-ad67-41aa-a965-a4f53ae68fb9" />

The reservoir state matrix in the frequency domain (Eq. 17) will be computed after taking the Fourier transformation of Eq. 14 and equating it to the Fourier transformation output signal. (To see details, check the paper)
<img width="1026" height="265" alt="image" src="https://github.com/user-attachments/assets/70249393-2128-4921-acb0-83e0801ed93c" />

The following optimization formulation based on Eq. 17 is developed to determine the optimal eigenvalues and output weights. 
<img width="575" height="550" alt="image" src="https://github.com/user-attachments/assets/3b6f0abb-caab-4253-8abe-ac2987f58a40" />




# Optimization Theory
The optimization problem is formulated as:
