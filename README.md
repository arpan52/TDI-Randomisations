# TDI-Randomisations

This repository contains the R source code and a simulation protocol document for all simulation studies conducted in the paper:

**_Non-Orthogonal Factorial Designs for Trials of Therapist-Delivered Interventions: Allowing for Therapist Capacity, Availability and Turnover_**

by **_Arpan Singh_**, **_Steven G. Gilmour_**, and **_Rebecca E. A. Walwyn_**.

The paper provides the methodological foundation for randomisation strategies used in trials of practitioner-delivered interventions (e.g., psychotherapy, physiotherapy, and surgery). It proposes three novel randomisation methods that explicitly account for practical constraints faced by practitioners (hereafter referred to as **therapists**), including therapist capacity, availability, and turnover.

As the number of simulations or the sample size increases, the computational burden also increases. Therefore, the use of a **High-Performance Computing (HPC)** facility is recommended.

These methods are compared with existing randomisation strategies in terms of:

- **Operational efficiency**, measured by the **time to target recruitment**, and
- **Statistical efficiency**, evaluated through the **precision of inference** obtained from the subsequent analysis.

All scripts required to reproduce the simulation results and figures presented in the paper are included in this repository.

**Paper:** LINK

---

# 📘 Description of the Code

All simulation studies in this repository are based on six randomisation strategies used to generate a trial design and two statistical models used to generate synthetic datasets from each of these designs. The generated data are then analysed to assess the performance of all six novel and existing randomisation methods. Details of these statistical models and their parameters can be found in the original paper.

---

# 🔧 Inputs

To run a simulation, the user must provide the following input parameters.

## Parameters for Randomisation Methods and Data Simulation

- `N`: Integer (>0) — Trial sample size
- `n_C`: Integer (>=1) — Number of centres
- `n_B`: Integer (>=1) — Number of batches
- `n_I`: Integer (>=1) — Number of interventions
- `E_n_T`: Integer (>=1) — Expected number of therapists per centre in a multi-centre trial
- `n_T_single`: Integer (>=1) — Number of therapists in a single-centre trial
- `Ex_C_P`: Numeric scalar (>=0) — Expected capacity of a part-time (embedded) therapist over the course of the trial
- `Ex_C_F`: Numeric scalar (>=0) — Expected capacity of a full-time (research) therapist over the course of the trial
- `sd_P`: Numeric scalar (>0) — Standard deviation of the capacity of a part-time therapist
- `sd_F`: Numeric scalar (>0) — Standard deviation of the capacity of a full-time therapist
- `a_P`: Numeric scalar between 0 and 1 — Proportion of part-time therapists in a centre
- `a`: Numeric scalar between 0 and 1 — Lower bound of the Uniform(`a`, 1) distribution used to generate the proportion of available therapists
- `frac`: Numeric scalar — Proportion of additional recruitment permitted within a batch under the flexible randomisation method
- `Add_Batch`: Binary (0/1) — Indicator of whether additional batches are added in the event of under-recruitment (1 = Yes, 0 = No)
- `n_sim`: Integer (>=1) — Number of simulation runs
- `time`: Binary (0/1) — Data-generating model (0 = Model 1, 1 = Model 2)

## Parameters for Data Analysis

### (1) Fixed-effect Parameters

- `del_0`: Numeric scalar — Overall mean
- `del_1`: Numeric scalar — Intervention effect
- `del_2`: Numeric scalar — Time effect (not applicable when `time = 0`)
- `del_3`: Numeric scalar — Time-by-intervention interaction effect (not applicable when `time = 0`)
- `del_4`: Numeric scalar — Patient-level matching covariate effect

### (2) Random-effect Parameters

- `sig2_T`: Numeric `n_C × 1` vector (>0) — Variance of therapist random effects for each centre
- `sig2_B`: Numeric scalar (>0) — Variance of the batch random effect (not applicable when `time = 1`)
- `sig2_C`: Numeric scalar (>0) — Variance of the centre random effect
- `sig2_IT`: Numeric scalar (>0) — Variance of the therapist-by-intervention random effect
- `sig2_IB`: Numeric scalar (>0) — Variance of the batch-by-intervention random effect (not applicable when `time = 1`)
- `sig2_IC`: Numeric scalar (>0) — Variance of the centre-by-intervention random effect
- `sig2_TB`: Numeric scalar (>0) — Variance of the therapist-by-batch random effect (not applicable when `time = 1`)
- `sig2_CB`: Numeric scalar (>0) — Variance of the centre-by-batch random effect (not applicable when `time = 1`)
- `sig2_ITB`: Numeric scalar (>0) — Variance of the therapist-by-batch-by-intervention random effect (not applicable when `time = 1`)
- `sig2_ICB`: Numeric scalar (>0) — Variance of the centre-by-batch-by-intervention random effect (not applicable when `time = 1`)
- `sig2_T_learn`: Numeric scalar (>0) — Variance of the therapist learning random effect (not applicable when `time = 0`)
- `sig2_C_learn`: Numeric scalar (>0) — Variance of the centre learning random effect (not applicable when `time = 0`)
- `sig2_IT_learn`: Numeric scalar (>0) — Variance of the therapist-by-intervention learning random effect (not applicable when `time = 0`)
- `sig2_IC_learn`: Numeric scalar (>0) — Variance of the centre-by-intervention learning random effect (not applicable when `time = 0`)
- `sig2_E`: Numeric scalar (>0) — Variance of the patient-specific residual error

---

# 📤 Reading the Output from *HPC_code_model_1_6_rands_multicentre_simulate_D_and_A.R*


The script produces the following two outputs.

## (1) `time_to_target_avg` (1 × 6 vector)

This vector contains the average values of the matrix `time_to_target`, which has dimensions **n_sim × 6**.

- Each **row** corresponds to one simulation run.
- Each **column** corresponds to one of the six randomisation methods described in the paper.

The *(i, j)*-th entry of `time_to_target` represents the number of recruitment batches required by the *j*-th randomisation method in the *i*-th simulation to achieve the target sample size `N`, while satisfying the therapist capacity, availability, and turnover constraints specified by the input parameters.

This quantity is referred to as the **Recruitment Time to Target** in the paper.

The reported output, `time_to_target_avg`, is obtained by averaging each column of `time_to_target` over the `n_sim` simulation runs.

A typical output is shown below.

| Non.rand | Wal25 | Flexible | SSPP2 | IBD | Flex.SPP2 |
|----------:|------:|---------:|------:|----:|----------:|
| 7.2 | 6.0 | 7.2 | 7.2 | 11.1 | 7.3 |

For example, the average recruitment time to target under the **IBD** randomisation strategy is **11.1** batches for the specified input parameter values.

---

### (2) `combined_standard` (6 × 16 matrix)

This matrix summarises the simulation results by reporting the average estimated model parameters together with several performance measures, averaged over the `n_sim` simulation runs for each randomisation strategy. For presentation purposes, the table below includes an additional row containing the true parameter values used to generate the simulated data.

| Randomisation | m | s | th | it | t2_er | bdry | fail | mse_marg | b | c | tb | ib  | cb | ic | itb | icb |
|:--------------|--:|--:|--:|--:|------:|-----:|-----:|---------:|--:|--:|--:|--:|--:|--:|--:|--:|
| **True Values** | 0.4000000 | 0.0000000 | 0.10000000 | 0.1500000 | 0.0000000 | 0.0000000 | 0.0 | 0.0000000 | 0.10000000 | 0 | 0.10000000 | 0.10000000 | 0 | 0 | 0.10000000 | 0 |
| **Non.rand** | 0.4442207 | 0.1413485 | 0.86743851 | 0.2130246 | 0.4000000 | 0.7000000 | 0.0 | 1.6207921 | 0.06890130 | NaN | 0.08830071 | 0.08504451 |NaN | NaN | 0.16093465 | NaN |
| **Wal25** | 0.3074686 | 0.1277925 | 0.11589185 | 0.1488394 | 0.4000000 | 0.3000000 | 0.0 | 0.9629645 | 0.07598271 | NaN | 0.09649396 | 0.10201022 |NaN | NaN | 0.12671330 | NaN |
| **Flexible** | 0.3418612 | 0.1343041 | 0.11480754 | 0.1574232 | 0.2000000 | 0.7000000 | 0.0 | 0.9700090 | 0.03433533 | NaN | 0.09380883 | 0.14164963 | NaN | NaN | 0.12847781 | NaN |
| **SPP2** | 0.3742867 | 0.1221458 | 0.12025883 | 0.1860772 | 0.3333333 | 0.2222222 | 0.1 | 0.9807015 | 0.06666675 | NaN | 0.10077534 | 0.08815796 |NaN | NaN | 0.09558719 | NaN |
| **IBD** | 0.4174948 | 0.1175016 | 0.05837238 | 0.1939579 | 0.1111111 | 0.5555556 | 0.1 | 0.9890080 | 0.12566384 | NaN | 0.08027903 | 0.08201259 | NaN | NaN | 0.13593354 | NaN |
| **Flex.SPP2** | 0.3702524 | 0.1282595 | 0.05369462 | 0.1914197 | 0.3000000 | 0.6000000 | 0.0 | 0.9561290 | 0.11008239 | NaN | 0.08507082 | 0.10367278 |NaN | NaN | 0.07517699 | NaN |


The output columns are described below.

| Column | Description |
|:------|:------------|
| `m` | Average estimate of `del_1` |
| `s` | Standard deviation of the estimate of `del_1` over the `n_sim` simulation runs |
| `th` | Average estimate of `sig2_T` |
| `it` | Average estimate of `sig2_IT` |
| `t2_er` | Type II error rate |
| `bdry` | Proportion of simulations in which one or more random-effect variance estimates were equal to zero when the corresponding true variance was non-zero |
| `fail` | Proportion of REML fitting failures |
| `mse_marg` | Average marginal mean squared error |
| `b` | Average estimate of `sig2_B` |
| `c` | Average estimate of `sig2_C` |
| `tb` | Average estimate of `sig2_TB` |
| `ib` | Average estimate of `sig2_IB` |
| `cb` | Average estimate of `sig2_CB` |
| `ic` | Average estimate of `sig2_IC` |
| `itb` | Average estimate of `sig2_ITB` |
| `icb` | Average estimate of `sig2_ICB` |

For example, under the **Non.rand** randomisation strategy, the average estimated therapist variance (`sig2_T`) is **0.8674**, whereas the true value is **0.10**, indicating substantial upward bias in the estimated therapist variance under this simulation scenario.
