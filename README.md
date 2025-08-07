# Sensory_Memory_Computational_Model
A computational model predicting changes to neural activity (adaptation) and to behavior (habituation) following repeated stimulation. 

The model is based in MATLAB (R2023a tested). Scripts are organized as follows and as depicted inthe schematic below:
* _Optimization.m_ generates optimized model parameters using the 20s/60s experimental AWA neural and behavioral response dataset (in _ExpData_conc.mat_) and the optimization cost function (_cost_fn.m_).
* _simulate_neural_activity_fn.m_ produces a simulated neural response timecourse based on the input model and stimulation parameters.
* _simulate_behavioral_activity_fn.m_ produces a simulated behavior probability timecourse based on the simulated neural response output (specifically, the hiddent latent representation).
* _summarize_behavior.m_ calculates adaptation and habituation indices.
<img width="1103" height="653" alt="image" src="https://github.com/user-attachments/assets/4a357d0d-e943-472a-80a5-cbb85dc41e37" />


To generate a new model prediction of neural acrtivity and response, see demo code _Demo_model_prediction.m_.
<img width="726" height="447" alt="image" src="https://github.com/user-attachments/assets/58a05e44-fded-4db8-b595-713f6b019bf6" />
<img width="726" height="428" alt="image" src="https://github.com/user-attachments/assets/c8dfc57e-c7dc-4710-946a-5fb3c5cf0d47" />
