# Fast_Decoupled_Power_Flow 9 Bus Case 

Fast Decoupled Power Flow (FDPF) is a numerical technique used in power system analysis to determine the steady-state operating conditions of an electrical network. It is a simplified version of the Gauss-Seidel method and is widely employed in power engineering for solving load flow problems. Load flow analysis is crucial in power system planning, operation, and control as it helps in understanding how power flows through the network under different operating conditions.

In traditional power flow methods, all bus variables (voltage magnitudes and angles) are updated simultaneously in each iteration. However, FDPF takes advantage of the inherent characteristics of power systems to decouple the equations governing voltage angles and magnitudes, allowing for a more efficient solution process. By decoupling these equations, FDPF significantly reduces the number of iterations required to reach a solution, making it computationally faster than the Gauss-Seidel method.

The decoupling of equations is based on the assumptions that the power system can be divided into two subsystems: the P-V buses (P: active power, V: voltage magnitude) and the Q-V buses (Q: reactive power, V: voltage magnitude). These assumptions simplify the complex interactions between different buses and allow for a more straightforward solution process.

In summary, Fast Decoupled Power Flow is a powerful algorithm used by power engineers to quickly and efficiently analyze and optimize power system operating conditions. By simplifying the iterative process, FDPF provides accurate results in a fraction of the time compared to traditional methods, making it an essential tool in the field of power system analysis and planning.

In Fast Decoupled Power Flow unlike the NR we are not consuming time for calculation of Jacobi, FDPF is fast and making more iterations in less time and it is making approximations instead of exact result.
>|Gij|=0 & |Vi|=0 & Sinθij=0 & Cosθ ij=1

## Results
![image](https://github.com/yavuzCodiin/Fast_Decoupled_Power-_Flow/assets/82445309/fcbd8a8f-fec7-46a9-9199-0cd6b368d152)
