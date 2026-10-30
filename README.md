# Brain-State-Resolved Consistency of Corticospinal Responses to TMS

MEP amplitudes vary substantially from trial to trial. This project
examines whether ongoing pre-stimulus EEG activity is associated not
only with MEP amplitude, but also with the consistency of corticospinal
responses across trials.

The analysis pipeline supports the extraction of pre-stimulation EEG
features, assignment of trials to brain-state-defined subsets, and
quantification of MEP variability within those subsets. The published
study investigated pre-stimulation sensorimotor mu dynamics, including
mu power, phase, and interhemispheric M1--M1 functional connectivity.

The analyses accompany:

> van Hattem, T., Hougland, J. R., Ahola, O., Goetz, S. M., Humaidan,
> D., Jooß, A., & Ziemann, U. (2026). Brain-state-resolved consistency
> of corticospinal responses to TMS. *Brain Stimulation, 19*(6), 103239.
> https://doi.org/10.1016/j.brs.2026.103239

## Running the analyses

Run the scripts in the order appropriate for your dataset and analysis:

1.  **Preprocess EEG and EMG** 
2.  **Extract MEP amplitudes**
3.  **Extract pre-stimulus EEG features**
4.  **Split trials by pre-stimulation features**
5.  **Quantify MEP variability**

The scripts contain user-adjustable settings for selecting the brain
region, EEG feature, frequency bands, trial-splitting method, and
output/plotting options. Review these settings before running an
analysis.

## Citation

If you use this code, please cite the associated publication:

van Hattem, T., Hougland, J. R., Ahola, O., Goetz, S. M., Humaidan, D.,
Jooß, A., & Ziemann, U. (2026). Brain-state-resolved consistency of
corticospinal responses to TMS. *Brain Stimulation, 19*(6), 103239.
https://doi.org/10.1016/j.brs.2026.103239

## Contact

For questions about the analysis code, please contact:

**Timo van Hattem**\
Brain Networks & Plasticity Lab\
Hertie Institute for Clinical Brain Research\
University of Tübingen\
timo.van-hattem@uni-tuebingen.de
