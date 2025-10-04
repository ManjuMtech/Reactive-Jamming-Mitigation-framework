# Cooperative Framework against Reactive Jamming in Wireless Systems

**Research Project (JTD792) | Jan 2025 – May 2025**  
**Under the guidance of:** *Dr. Harshan Jagadeesh*  
**Grade:** A+
---
## 🧠 Overview
This project investigates **cooperative communication techniques** to mitigate **reactive jamming attacks** in wireless systems.  
The framework simulates a **two-user cooperative setup** involving:
- A **victim transmitter** using On–Off Keying (OOK)
- A **helper node** using M-QAM modulation
- A **joint Maximum A Posteriori (MAP) decoder** at the receiver.

The system is evaluated under varying **energy division factors (α)**, **SNR values**, and **number of receiver antennas** to study performance degradation and cooperative benefits.
---
## ⚙️ Simulation Details
### System Components
| Component | Description |
|------------|-------------|
| Victim Transmitter | OOK-modulated signal |
| Helper Node | M-QAM (e.g., 16-QAM) cooperative transmission |
| Receiver | Joint MAP decoder combining both signal paths |
| Jammer Model | Reactive jammer activating upon energy detection |
| Channel Model | AWGN with optional fading components |

### Parameters
- **Energy division factor (α):** Controls power sharing between victim and helper  
- **SNR range:** 20 dB – 35 dB  
- **Helper receive antennas:** 1 to 4  
- **Performance metric:** Bit Error Rate (BER)
---
## 📊 Results
- **BER vs. α Curves:**  
  For a given SNR, BER increases with higher α — indicating that allocating more energy to the victim link makes the system more vulnerable to jamming.  
- **Effect of SNR and antennas:**  Pe_vs_alpha @Nc_2.png
  Increasing SNR or adding more helper antennas reduces BER significantly, improving robustness under jamming conditions.

Example plot:

