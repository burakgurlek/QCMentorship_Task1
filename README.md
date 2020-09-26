# QCMentorship_Task1

Has been created for solution of Task 1 in QC Mentorship program's selection process. This has been implemented by using Cirq (0.8.2). Since Cirq is currently in the development please use the referred version above.

Mostly, the minumum distance is found for a circuit with one layer. In the case of multiple layers, my algorith is not convergent. I guess the reason is wrong selection of the optimization method or its parameters. One needs to be careful since gradient of the distance is not smooth.
To resolve this, one can use grid search but I did not implemented this here due to time restrictions. I believe circuit itself is correctly implemented. 

