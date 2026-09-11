% Description: return the ADC sample mode type for Verasonics hardware.
% Input is the sample mode number: 50, 67, 100, 200%

function ADC_sampleMode = getADCSampleMode(ADC_sampleModeNumber)
    switch ADC_sampleModeNumber
        case 50
            ADC_sampleMode = 'BS50BW';
        case 67
            ADC_sampleMode = 'BS67BW';
        case 100
            ADC_sampleMode = 'BS100BW';
        case 200
            ADC_sampleMode = 'NS200BW';
    end