# Chapter 1: Fundamentals of Vibration

Text: Fundamentals of Acoustics
Authors: Lawrence E. Kinsler, Austin R. Frey, Alan B. Coppens, and James V. Sanders
Year Published: 2000

## 1.1 Introduction

The chapter begins by establishing the need for a consistent system of units in acoustics. The field encompasses a wide range of scientific and engineering disciplines, leading to a lack of uniformity in unit systems across the literature. Early work used the CGS (centimeter-gram-second) system, while engineering applications often employed a mixture of metric and English units. Electroacoustics and underwater acoustics commonly used the MKS (meter-kilogram-second) system.

The text adopts the SI (Le Système International d'Unités) system as the standard, which is a codification of the MKS system. Throughout the text, "log" represents logarithm to base 10, while "ln" represents the natural logarithm (logarithm to base e).

*Acoustics is defined as the science of generation, transmission, and reception of energy as vibrational waves in matter.* When molecules of a fluid or solid are displaced from their normal configurations, an internal elastic restoring force arises. This elastic restoring force, coupled with the inertia of the system, enables matter to participate in oscillatory vibrations and thereby generate and transmit acoustic waves.

Examples of elastic restoring forces include:

- Tensile force when a spring is stretched
- Increased pressure when a fluid is compressed  
- Restoring force when a point on a stretched wire is displaced transverse to its length

The most familiar acoustic phenomenon is sound, perceived when vibrational disturbances have frequencies in the range from about 20 Hz to 20,000 Hz (where 1 Hz = 1 cycle per second). However, acoustics in a broader sense includes ultrasonic frequencies above 20,000 Hz and infrasonic frequencies below 20 Hz. The nature of vibrations associated with acoustics encompasses many types, from simple sinusoidal vibrations of a tuning fork to complex vibrations of a bowed violin string and nonperiodic motions in explosions.

/*

This introduction provides a few insights that mark out boundaries for a study of acoustics oriented toward digital sound engineering and audio synthesis. 

First, its simple definition of acoustics: *the science of the generation, transmission, and reception of energy as vibrational waves in matter.* 

And second, its statement of the range of frequencies over which acoustic vibrations can be detected as sound: 20 Hz to 20,000 Hz.

These axioms lead to a few key questions:

1. How do we percieve sound between 20 Hz and 20,000 Hz? I am familiar with the scales, keys, and notes of classical Western music -- so this question may be easiest to answer through a proxy: By what patterns do sounds corresponding to these frequencies match to the pitches, scales, chords of classical Western music? 

2. What is the relationship between binary descriptions of sound and analog acoustic phenomena? How is binary audio data (digital descriptions of sound waves) coupled to analog sound output devices? I imagine the fundamental relationship is not dissimilar from the relationship between binary (or digitized) descriptions of optics and 'continuous' classical optical phenomena... Which, I suppose, leads to another question. 

3. Three main theoretical frameworks are used to describe light: 

        a. Geometric optics, wherein light is described as sets of rays
        b. Physical optics, wherein light is described as continuous waves
        c. Quantum optics, wherein light is described as photons which exhibit both partical-like and wave-like characteristics

    What are the main theoretical frameworks used to describe acoustic phenomena?

*/

## 1.2 The Simple Oscillator

A mass $m$ fastened to a spring and constrained to move parallel to the spring represents the fundamental model. When displaced slightly from rest and released, the mass vibrates with displacement that is a sinusoidal function of time. These vibrations are called simple harmonic vibrations.

Many acoustic vibrators can be modeled as simple oscillators, including loaded tuning forks and loudspeaker diaphragms. More complex vibrating systems often exhibit characteristics of simple systems and can be approximated by simple oscillators as a first approach.

The physical restrictions for a simple oscillator are:
- The restoring force must be directly proportional to displacement (Hooke's law)
- The mass must be constant
- There must be no losses to attenuate the motion

Under these conditions, the frequency of vibration is independent of amplitude, and the motion is simple harmonic.

The restoring force $f$ in newtons (N) is expressed by:

$$f = -sx \quad (1.2.1)$$

where $x$ is displacement in meters (m) from rest position, $s$ is the stiffness or spring constant in N/m, and the minus sign indicates the force opposes the displacement.

Substituting into Newton's second law:

$$f = m\frac{d^2x}{dt^2} \quad (1.2.2)$$

yields:

$$\frac{d^2x}{dt^2} + \frac{s}{m}x = 0 \quad (1.2.3)$$

Defining the constant:

$$\omega_0^2 = s/m \quad (1.2.4)$$

the equation becomes:

$$\frac{d^2x}{dt^2} + \omega_0^2x = 0 \quad (1.2.5)$$

This important linear differential equation has a well-known general solution. Trial solutions of the form $x = A_1\cos\gamma t$ and $x = A_2\sin\omega_0 t$ can be verified by substitution. The complete general solution is:

$$x = A_1\cos\omega_0 t + A_2\sin\omega_0 t \quad (1.2.8)$$

where $A_1$ and $A_2$ are arbitrary constants, and $\omega_0$ is the natural angular frequency in radians per second (rad/s). The natural frequency $f_0$ in hertz (Hz) relates to the natural angular frequency by:

$$f_0 = \omega_0/2\pi \quad (1.2.9)$$

The period $T$ of one complete vibration is:

$$T = 1/f_0 \quad (1.2.10)$$

Note that decreasing stiffness or increasing mass lowers the frequency.

## 1.3 Initial Conditions

The arbitrary constants $A_1$ and $A_2$ are determined by initial conditions. If at time $t = 0$ the mass has initial displacement $x_0$ and initial speed $u_0$, direct substitution shows that $A_1$ equals the initial displacement $x_0$, and from differentiation of (1.2.8), $u_0 = \omega_0A_2$, so (1.2.8) becomes:

$$x = x_0\cos\omega_0 t + (u_0/\omega_0)\sin\omega_0 t \quad (1.3.1)$$

An alternative form can be obtained by letting $A_1 = A\cos\phi$ and $A_2 = -A\sin\phi$, where $A$ and $\phi$ are new arbitrary constants:

$$x = A\cos(\omega_0 t + \phi) \quad (1.3.2)$$

where $A$ is the amplitude of the motion and $\phi$ is the initial phase angle. The values are determined by initial conditions:

$$A = [x_0^2 + (u_0/\omega_0)^2]^{1/2} \quad \text{and} \quad \phi = \tan^{-1}(-u_0/\omega_0x_0) \quad (1.3.3)$$

Successive differentiation shows the speed is:

$$u = -U\sin(\omega_0 t + \phi) \quad (1.3.4)$$

where $U = \omega_0A$ is the speed amplitude, and the acceleration is:

$$a = -\omega_0U\cos(\omega_0 t + \phi) \quad (1.3.5)$$

The displacement lags 90° (π/2 rad) behind the speed, and the acceleration is 180° (π rad) out of phase with the displacement.

## 1.4 Energy of Vibration

The mechanical energy $E$ of the system is the sum of potential energy $E_p$ and kinetic energy $E_k$. The potential energy is the work done in distorting the spring:

$$E_p = \int_0^x sx\,dx = \frac{1}{2}sx^2 \quad (1.4.1)$$

Using the displacement expression (1.3.2):

$$E_p = \frac{1}{2}sA^2\cos^2(\omega_0 t + \phi) \quad (1.4.2)$$

The kinetic energy is:

$$E_k = \frac{1}{2}mu^2 \quad (1.4.3)$$

Using the speed expression (1.3.4):

$$E_k = \frac{1}{2}mU^2\sin^2(\omega_0 t + \phi) \quad (1.4.4)$$

The total energy is:

$$E = E_p + E_k = \frac{1}{2}m\omega_0^2A^2 \quad (1.4.5)$$

using $s = m\omega_0^2$, $U = \omega_0A$, and the identity $\sin^2\sigma + \cos^2\sigma = 1$. This can also be written as:

$$E = \frac{1}{2}sA^2 = \frac{1}{2}mU^2 \quad (1.4.6)$$

The total energy is constant (independent of time) and equals either the maximum potential energy (when the mass is at greatest displacement and instantaneously at rest) or the maximum kinetic energy (when the mass passes through equilibrium with maximum speed). Since no external forces or friction were assumed, energy conservation is satisfied.

## 1.5 Complex Exponential Method of Solution

The text adopts the engineering convention of representing time dependence of oscillatory functions by $\exp(j\omega t)$ rather than the physics convention of $\exp(-i\omega t)$, due to analogies with engineering applications. This may sometimes result in exchange of complex functions, but context usually resolves ambiguities.

A more general approach to solving linear differential equations like (1.2.5) is to postulate:

$$x = Ae^{\gamma t} \quad (1.5.1)$$

Substitution gives $\gamma^2 = -\omega_0^2$ or $\gamma = \pm j\omega_0$. Thus, the general solution is:

$$x = A_1e^{j\omega_0 t} + A_2e^{-j\omega_0 t} \quad (1.5.2)$$

where $A_1$ and $A_2$ are determined by initial conditions $x(0) = x_0$ and $dx(0)/dt = u_0$:

$$A_1 + A_2 = x_0 \quad \text{and} \quad A_1 - A_2 = u_0/j\omega_0 = -ju_0/\omega_0 \quad (1.5.3)$$

Solving:

$$A_1 = \frac{1}{2}(x_0 - ju_0/\omega_0) \quad \text{and} \quad A_2 = \frac{1}{2}(x_0 + ju_0/\omega_0) \quad (1.5.4)$$

Since $A_1$ and $A_2$ are complex conjugates, the real solution is:

$$x = x_0\cos\omega_0 t + (u_0/\omega_0)\sin\omega_0 t \quad (1.5.5)$$

which is identical to (1.3.1).

For the complex solution, if we express $A_1 = a_1 + jb_1$ and $A_2 = a_2 + jb_2$ in (1.5.2), taking the real part gives:

$$\text{Re}\{x\} = (a_1 + a_2)\cos\omega_0 t - (b_1 - b_2)\sin\omega_0 t \quad (1.5.6)$$

A complete solution is obtained if displacement is written in complex form:

$$x = Ae^{j\omega_0 t} \quad (1.5.7)$$

where $A = a + jb$, and only the real part is considered:

$$\text{Re}\{x\} = a\cos\omega_0 t - b\sin\omega_0 t \quad (1.5.8)$$

From form (1.5.7), the complex speed is $u = j\omega_0Ae^{j\omega_0 t} = j\omega_0x$ and complex acceleration is $a = -\omega_0^2Ae^{j\omega_0 t} = -\omega_0^2x$.

The expression $\exp(j\omega_0 t)$ can be thought of as a phasor of unit length rotating counterclockwise in the complex plane with angular speed $\omega_0$. Any complex quantity $A = a + jb$ can be represented by a phasor of length $A = \sqrt{a^2 + b^2}$ making angle $\phi = \tan^{-1}(b/a)$ counterclockwise from the positive real axis. The product $A\exp(j\omega_0 t)$ represents a phasor of length $A$ and initial phase angle $\phi$ rotating with angular speed $\omega_0$. The real part of this rotating phasor (its projection on the real axis) is:

$$A\cos(\omega_0 t + \phi) \quad (1.5.11)$$

and varies harmonically with time.

## 1.6 Damped Oscillations

When a real body oscillates, dissipative (frictional) forces arise, resulting in damping of the oscillations—a decrease in amplitude of free oscillations with time. A viscous frictional force $f_r$ on a simple oscillator, assumed proportional to speed and opposing motion, can be expressed as:

$$f_r = -R_m\frac{dx}{dt} \quad (1.6.1)$$

where $R_m$ is a positive constant called the mechanical resistance with units of newton-second per meter (N⋅s/m) or kilogram per second (kg/s).

Including resistance in the equation of motion:

$$m\frac{d^2x}{dt^2} + R_m\frac{dx}{dt} + sx = 0 \quad (1.6.2)$$

Dividing by $m$ and recalling $\omega_0 = \sqrt{s/m}$:

$$\frac{d^2x}{dt^2} + \frac{R_m}{m}\frac{dx}{dt} + \omega_0^2x = 0 \quad (1.6.3)$$

Using the complex exponential method with $x = Ae^{\gamma t}$:

$$\gamma^2 + (R_m/m)\gamma + \omega_0^2 = 0 \quad (1.6.5)$$

This gives:

$$\gamma = -\beta \pm (\beta^2 - \omega_0^2)^{1/2} \quad (1.6.7)$$

where:

$$\beta = R_m/2m \quad (1.6.8)$$

In most cases of acoustic importance, the mechanical resistance $R_m$ is small enough that $\omega_0 > \beta$ and $\gamma$ is complex. Defining a new constant:

$$\omega_d = (\omega_0^2 - \beta^2)^{1/2} \quad (1.6.10)$$

Then:

$$\gamma = -\beta \pm j\omega_d \quad (1.6.11)$$

where $\omega_d$ is the natural angular frequency of the damped oscillator. Note that $\omega_d$ is always less than $\omega_0$ for the same oscillator without damping.

The complete solution is:

$$x = e^{-\beta t}(A_1e^{j\omega_d t} + A_2e^{-j\omega_d t}) \quad (1.6.12)$$

One convenient form of the general solution is:

$$x = Ae^{-\beta t}\cos(\omega_d t + \phi) \quad (1.6.13)$$

where $A$ and $\phi$ are real constants determined by initial conditions.

The amplitude of the damped oscillator, defined as $A\exp(-\beta t)$, is no longer constant but decreases exponentially with time. The frequency is independent of oscillation amplitude.

The time required for amplitude to decrease to $1/e$ of its initial value is the relaxation time (also called decay modulus, decay time, time constant, or characteristic time):

$$\tau = 1/\beta = 2m/R_m \quad (1.6.14)$$

The quantity $\beta$ is the temporal absorption coefficient. The smaller $R_m$, the larger $\tau$ and the longer it takes for oscillations to damp out.

If mechanical resistance $R_m$ is large enough that $\omega_0 \leq \beta$, the system is no longer oscillatory; a displaced mass returns asymptotically to rest. If $\beta = \omega_0$, the system is critically damped.
