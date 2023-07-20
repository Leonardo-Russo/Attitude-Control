# Attitude Control Study

This repository is dedicated to the study and implementation of different methods for attitude control. The projects focus on three prominent methodologies: Quasi-Rigid Body Theory, Kane's Method, and Non-linear Attitude Control.

## Projects

### 1. Quasi-Rigid Body Theory

Quasi-rigid body theory is a simplification of the motion dynamics of systems that deform slightly under loads, which can be effectively treated as a rigid body for certain analyses.

#### Theory

In Quasi-Rigid Body Theory, the complex dynamics of deformable bodies are approximated by rigid body dynamics, plus a small number of additional generalized coordinates to represent the deformation. It provides a balance between modeling fidelity and computational efficiency, allowing more complex dynamics to be captured than rigid body dynamics alone, while remaining tractable for analysis and control design.

### 2. Kane's Method

Kane's method is a powerful tool in dynamics, and it's used to derive equations of motion for complex mechanical systems.

#### Theory

Kane's method employs generalized speeds and non-redundant forces to derive the equations of motion. It avoids the need to eliminate constraint forces from the equations, simplifying the process considerably for complex systems. In the context of attitude control, it can be used to derive the equations of motion for a spacecraft and design a control system to maintain or change its attitude.

### 3. Non-linear Attitude Control

Non-linear control techniques can account for the non-linearities inherent in the attitude dynamics of a spacecraft and provide more accurate and robust control than linear methods.

#### Theory

Non-linear attitude control involves designing a control law that drives the current state of the system towards a desired state. This usually involves the computation of a desired torque based on the current and desired attitudes and angular velocities. It can account for the non-linearities in the attitude dynamics and constraints on the control input, leading to more accurate and robust control performance.

---

## Disclaimer

These projects are created for educational purposes. Please use responsibly and ensure all simulations and data are used in a manner adhering to applicable laws and regulations.

