# deim-cur-ecg
This repository contains research code written in MATLAB for applying DEIM CUR with incremental QR to synthetic and real ECG data.
The code corresponds to the analyses and results described in a submitted paper:
“Finding Representative Electrocardiogram Beat Morphologies with CUR.” Hendryx, E., Rusin, C., Riviere, B., Sorensen, D. 
Submitted 2016.

The generation and analysis of the synthetic data can be carried out without dependencies on other code sources. In addition, the specific synthetic data sets constructed for controled sensitivity testing in the above mentioned paper is stored here. 

The real data studied in this repository includes the MIT-BIH Arrhythmia Database [1], the Massachusets General Hospital-Marquette Foundation (MGH-MF) Waveform Database [2], and the St. Petersburg Institute of Cardiological Technics (Incart) 12-lead Arrhythmia Database. All three of these data sets are availble through PhysioNet [3]. In addition to needing the .mat, .info, and rr.txt files from each record of these databases, the matrix constructor codes require the user to incorpate code from the plotATM.m function [4]. Because the MGH-MF records are longer, it is recommended that users use the wfdb2mat.m function [5] in WFDB Toolbox for Matlab and Octave [6] to download the corresponding .mat files. Both plotATM.m and the WFDB Toolbox for Matlab and Octave are avaliable on PhysioNet [3].

References:

[1] Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)

[2] Welch JP, Ford PJ, Teplick RS, Rubsamen RM. The Massachusetts General Hospital-Marquette Foundation Hemodynamic and Electrocardiographic Database -- Comprehensive collection of critical care waveforms. J Clinical Monitoring 7(1):96-97 (1991).

[3] Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages; http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13).

[4] Abdala, O., Hislop, H. (2014). "plotATM.m" (Version 1.1) [Computer code].	Available at https://physionet.org/cgi-bin/atm/ATM (Accessed March 30, 2016.)

[5] Silva, I. (2014). "wfdb2mat.m" (Version 0.1) [Computer code]. Available at https://physionet.org/physiotools/matlab/wfdb-app-matlab/ (Accessed May 30, 2017).

[6] Silva, I, Moody, G. "An Open-source Toolbox for Analysing and Processing PhysioNet Databases in MATLAB and Octave." Journal of Open Research Software 2(1):e27 [http://dx.doi.org/10.5334/jors.bi] ; 2014 (September 24).
