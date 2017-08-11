# deim-cur-ecg
This repository contains research code written in MATLAB for applying DEIM CUR with incremental QR to synthetic and real ECG data.
The code corresponds to the analyses and results described in a submitted paper:
“Finding Representative Electrocardiogram Beat Morphologies with CUR.” Hendryx, E., Rusin, C., Riviere, B., Sorensen, D. 
Submitted 2016.

The generation and analysis of the synthetic data can be carried out without dependencies on other code sources. In addition, the specific synthetic data sets constructed for controled sensitivity testing in the above mentioned paper is stored here. 

The real data studied in this repository includes the MIT-BIH Arrhythmia Database [1], the Massachusets General Hospital-Marquette Foundation (MGH-MF) Waveform Database [2], the St. Petersburg Institute of Cardiological Technics (Incart) 12-lead Arrhythmia Database, and the MIT-BIH Noise Stress Test Database (NSTDB) [3]. All three of these data sets are availble through PhysioNet [4]. In addition to needing the .mat, .info, and rr.txt files from each record of these databases, the matrix constructor codes require the user to incorpate code from the plotATM.m function [5]. Because the MGH-MF records are longer, it is recommended that users use the wfdb2mat.m function [6] in WFDB Toolbox for Matlab and Octave [7] to download the corresponding .mat files. Both plotATM.m and the WFDB Toolbox for Matlab and Octave are avaliable on PhysioNet [4].

While there is separate code included for each of these data sets, some analyses on the DS1 and DS2 subsets of the MIT-BIH Arrhythmia database (without added noise) and the first 40 files only of the MGH-MF Waveform Database are not currently included. These analyses can, however, be created by editing related codes in this repository to read in only the corresponding subset of the data files currently analyzed. Note that some scripts load previously saved result files, so users should carefully track the names of files to which results are saved as well as when these files are used for further analysis.


References:

[1] Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)

[2] Welch JP, Ford PJ, Teplick RS, Rubsamen RM. The Massachusetts General Hospital-Marquette Foundation Hemodynamic and Electrocardiographic Database -- Comprehensive collection of critical care waveforms. J Clinical Monitoring 7(1):96-97 (1991).

[3] G. B. Moody, W. Muldrow, R. G. Mark, A noise stress test for arrhythmia detectors, Computers in cardiology 11 (3) (1984) 381–384.

[4] Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages; http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13).

[5] Abdala, O., Hislop, H. (2014). "plotATM.m" (Version 1.1) [Computer code].	Available at https://physionet.org/cgi-bin/atm/ATM (Accessed March 30, 2016.)

[6] Silva, I. (2014). "wfdb2mat.m" (Version 0.1) [Computer code]. Available at https://physionet.org/physiotools/matlab/wfdb-app-matlab/ (Accessed May 30, 2017).

[7] Silva, I, Moody, G. "An Open-source Toolbox for Analysing and Processing PhysioNet Databases in MATLAB and Octave." Journal of Open Research Software 2(1):e27 [http://dx.doi.org/10.5334/jors.bi] ; 2014 (September 24).

