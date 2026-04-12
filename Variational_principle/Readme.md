# 📘 1D Periodic Potential — Plane Wave Method

This project solves the Schrödinger equation for an electron in a **1D periodic potential** using a **plane-wave variational approach**.

---

## 🔹 Problem Setup

We solve the time-independent Schrödinger equation:

$$
\left(-\frac{1}{2}\nabla^2 + V(x)\right)\psi(x) = E \psi(x)
$$

with a periodic potential:

$$
V(x + a) = V(x)
$$

We choose:

$$
V(x) = -V_0 \left[1 + \cos\left(\frac{2\pi x}{a}\right)\right]
$$

Using **Bloch’s theorem**, the wavefunction is written as:

$$
\psi_k(x) = \frac{1}{\sqrt{a}} \sum_G c_G \, e^{i(k+G)x}
$$

---

## 🔹 Goal

Convert the Schrödinger equation into a **matrix eigenvalue problem**:

$$
H \mathbf{c} = E \mathbf{c}
$$

---

## 🔹 Hamiltonian Matrix Elements (Derivation)

We compute matrix elements in the plane-wave basis:

$$
H_{G,G'} = \langle \phi_G | \hat{H} | \phi_{G'} \rangle
$$

where:

$$
\phi_G(x) = \frac{1}{\sqrt{a}} e^{iGx}
$$

---

## ✅ Step 1: Kinetic Energy Term

The kinetic operator is:

$$
\hat{T} = -\frac{1}{2}\nabla^2
$$

Acting on a plane wave:

$$
-\frac{1}{2}\nabla^2 e^{i(k+G)x}
= \frac{1}{2}(k+G)^2 e^{i(k+G)x}
$$

Thus, the matrix element becomes:

$$
T_{G,G'} = \frac{1}{2}(k+G)^2 \, \delta_{G,G'}
$$

✔ **Diagonal in plane-wave basis**

---

## ✅ Step 2: Potential Energy Term

We compute:

$$
V_{G,G'} = \frac{1}{a} \int_0^a V(x)\, e^{-i(G-G')x}\, dx
$$

Substitute the potential:

$$
V(x) = -V_0 - V_0 \cos\left(\frac{2\pi x}{a}\right)
$$

---

### 🔸 (i) Constant Term
  $$-V_0 \delta_{G,G'}$$

---

### 🔸 (ii) Cosine Term

Using:

$$
\cos\left(\frac{2\pi x}{a}\right)
= \frac{1}{2}\left(e^{i\frac{2\pi x}{a}} + e^{-i\frac{2\pi x}{a}}\right)
$$

we obtain coupling between plane waves:

$$
V_{G,G'} =
-\frac{V_0}{2}
\quad \text{if } G - G' = \pm \frac{2\pi}{a}
$$

---

## ✅ Final Hamiltonian

Combining both contributions:


$$H_{G,G'}(k) =
\frac{1}{2}(k+G)^2 \, \delta_{G,G'}- V_0 \, \delta_{G,G'}- \frac{V_0}{2} \, \delta_{G-G', \pm \frac{2\pi}{a}}$$


---

## 🔹 Key Insight

- **Diagonal terms** → kinetic energy + constant potential  
- **Off-diagonal terms** → coupling between plane waves  
- This coupling is responsible for:
  - Band formation  
  - Band gaps at Brillouin zone boundaries  
