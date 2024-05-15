# Conventions

VLBI and radio astronomy has many non-standard conventions when coming from physics. Additionally, these conventions change from telescope to telescope, often making it difficult 
to know what assumptions different data sets and codes are making. We will detail the 
specific conventions that Comrade adheres to.


## Rotation Convention

We follow the standard EHT and rotate starting from the upper y-axis and moving in a counter-clockwise direction. 

!!! note
    We still use the standard astronomy definition where the positive x-axis is to the left.


## Fourier Transform Convention

We use the positive exponent definition of the Fourier transform to define our visibilities. That is, we assume that the visibilities measured by a perfect interferometer are given by
```math
 V(u, v) = \int I(x, y)e^{2\pi i(ux + vy)}dx dy.
```
This convention is consistent with the AIPS convention and what is used in other EHT codes, such as eht-imaging. 

!!! warning
    This is the opposite convention of what is written in the EHT papers, but it is the correct version for the released data.


## Coherency matrix Convention

We use the factor of 2 definition when defining the coherency matrices. That is, the relation coherency matrix `C` is given by

```math
  C_{pq} = 
  2\begin{pmatrix}
    \left<v_{pa} v_{qa}^*\right> & \left< v_{pa}v_{qb}^*\right >\\
    \left<v_{pb} v_{qa}^*\right> & \left< v_{pb}v_{qb}^*\right >
  \end{pmatrix}.
```

where $v_{pa}$ is the voltage measured from sites $p$ and feed $a$.

#### Circular Polarization Conversions

To convert from measured $R/L$ circular cross-correlation products to the Fourier transform of the Stokes parameters, we use:

```math
  \begin{pmatrix}
      \tilde{I}\\ \tilde{Q} \\ \tilde{U} \\ \tilde{V}
  \end{pmatrix}
  =\frac{1}{2}
  \begin{pmatrix}
     \left<RR^*\right> + \left<LL^*\right> \\
     \left<RL^*\right> + \left<LR^*\right> \\
     i(\left<LR^*\right> - \left<RL^*\right>)\\
     \left<RR^*\right> - \left<LL^*\right>
  \end{pmatrix},
```

where e.g., $\left<RL^*\right> = 2\left<v_{pR}v^*_{pL}\right>$.

The inverse transformation is then:

```math
  C = 
  \begin{pmatrix}
     \tilde{I} + \tilde{V}  & \tilde{Q} + i\tilde{U}\\
     \tilde{Q} - i\tilde{U} & \tilde{I} - \tilde{V}
  \end{pmatrix}.
```

#### Linear Polarization Conversions

To convert from measured $X/Y$ linear cross-correlation products to the Fourier transform of the Stokes parameters, we use:

```math
  \begin{pmatrix}
      \tilde{I}\\ \tilde{Q} \\ \tilde{U} \\ \tilde{V}
  \end{pmatrix}
  =\frac{1}{2}
  \begin{pmatrix}
     \left<XX^*\right> + \left<YY^*\right> \\
     \left<XY^*\right> + \left<YX^*\right> \\
     i(\left<YX^*\right> - \left<XY^*\right>)\\
     \left<XX^*\right> - \left<YY^*\right>
  \end{pmatrix}.
```

The inverse transformation is then:

```math
  C = 
  \begin{pmatrix}
     \tilde{I} + \tilde{Q}  & \tilde{U} + i\tilde{V}\\
     \tilde{U} - i\tilde{V} & \tilde{I} - \tilde{Q}
  \end{pmatrix}.
```

where e.g., $\left<XY^*\right> = 2\left<v_{pX}v^*_{pY}\right>$.

