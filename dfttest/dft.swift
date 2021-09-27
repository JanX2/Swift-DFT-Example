//
//  dft.swift
//  dfttest
//
//  Created by Christopher Helf on 17.08.15.
//  Copyright (c) 2015-2021 Christopher Helf. All rights reserved.
//  Copyright (c) 2021-Present Jan Weiß. All rights reserved.

import Foundation
import Accelerate

class DFT {
    typealias Sample = Double
    typealias Complex = DSPDoubleComplex
    typealias SplitComplex = DSPDoubleSplitComplex
    
    fileprivate func getFrequencies(_ N: Int, samplerate: Sample) -> [Sample] {
        // Create an Array with the Frequencies
        let freqs = (0..<N/2).map {
            samplerate/Sample(N) * Sample($0)
        }
        
        return freqs
    }
    
    fileprivate func generateBandPassFilter(_ freqs: [Sample]) -> ([Sample], Int, Int) {
        var minIdx = freqs.count+1
        var maxIdx = -1
        
        let bandPassFilter: [Sample] = freqs.map {
            if ($0 >= self.lowerFreq && $0 <= self.higherFreq) {
                return 1.0
            } else {
                return 0.0
            }
        }
        
        for (i, element) in bandPassFilter.enumerated() {
            if (element == 1.0) {
                if(i<minIdx || minIdx == freqs.count+1) {
                    minIdx=i
                }
                if(i>maxIdx || maxIdx == -1) {
                    maxIdx=i
                }
            }
        }
        
        assert(maxIdx != -1)
        assert(minIdx != freqs.count+1)
        
        return (bandPassFilter, minIdx, maxIdx)
    }
    
    func calculate(_ inRealValues: [Sample], samplerate: Sample) {
        // ----------------------------------------------------------------
        // Size Variables
        // ----------------------------------------------------------------
        let N = inRealValues.count
        let halfN = N/2
        let N_ = vDSP_Length(N)
        let halfN_ = vDSP_Length(halfN)
        
        // ----------------------------------------------------------------
        // DFT & Variables Setup
        // ----------------------------------------------------------------
        let dftForwardSetup: vDSP_DFT_Setup = vDSP_DFT_zrop_CreateSetupD(nil, N_, .FORWARD)!
        let dftInverseSetup: vDSP_DFT_Setup = vDSP_DFT_zrop_CreateSetupD(dftForwardSetup, N_, .INVERSE)!
        
        // We need the complex buffer in two different data layouts!
        var tempComplex : [Complex] = [Complex](repeating: Complex(), count: halfN)
        
        var tempReal = [Sample](repeating: 0.0, count:halfN)
        var tempImag = [Sample](repeating: 0.0, count:halfN)
        
        // Split even and odd indexes of our input values into their own arrays
        // https://developer.apple.com/documentation/accelerate/performing_fourier_transforms_on_interleaved-complex_data
        tempReal.withUnsafeMutableBufferPointer { tempRealPtr in
            tempImag.withUnsafeMutableBufferPointer { tempImagPtr in
                
                var tempSplitComplex =
                    SplitComplex(realp: tempRealPtr.baseAddress!,
                                 imagp: tempImagPtr.baseAddress!)
                
                var valuesAsComplex: UnsafePointer<Complex>? = nil
                
                inRealValues.withUnsafeBytes {
                    valuesAsComplex =
                        $0.baseAddress?.bindMemory(to: Complex.self,
                                                   capacity: N)
                }
                
                vDSP_ctozD(valuesAsComplex!, vDSP_Stride(2),
                           &tempSplitComplex, vDSP_Stride(1),
                           halfN_)
            }
        }
        
        // ----------------------------------------------------------------
        // Forward DFT
        // ----------------------------------------------------------------
        
        // Do real->complex forward in-place DFT
        vDSP_DFT_ExecuteD(dftForwardSetup, tempReal, tempImag, &tempReal, &tempImag)
        
        // ----------------------------------------------------------------
        // Get the Frequency Spectrum
        // ----------------------------------------------------------------
        var fullSpectrum = [Sample](repeating: 0.0, count: halfN)
        
        // For polar coordinates
        var mag: [Sample] = [Sample](repeating: 0.0, count: halfN)
        var phase: [Sample] = [Sample](repeating: 0.0, count: halfN)
        
        tempReal.withUnsafeMutableBufferPointer { tempRealPtr in
            tempImag.withUnsafeMutableBufferPointer { tempImagPtr in
                
                var tempSplitComplex = SplitComplex(realp: tempRealPtr.baseAddress!,
                                                    imagp: tempImagPtr.baseAddress!)
                
                var dftMagnitudes = [Sample](repeating: 0.0, count: halfN)
                vDSP_zvmagsD(&tempSplitComplex, 1, &dftMagnitudes, 1, halfN_);
                
                // vDSP_zvmagsD returns squares of the DFT magnitudes, so take the root here
                let roots = sqrt(dftMagnitudes)
                
                // Normalize the Amplitudes
                vDSP_vsmulD(roots, vDSP_Stride(1), [1.0 / Sample(N)], &fullSpectrum, 1, halfN_)
                
                // ----------------------------------------------------------------
                // Convert from complex/rectangular (real, imaginary) coordinates
                // to polar (magnitude and phase) coordinates.
                // ----------------------------------------------------------------
                
                vDSP_zvabsD(&tempSplitComplex, vDSP_Stride(1),
                            &mag, vDSP_Stride(1),
                            halfN_);
                
                // Beware: phase values in output between -PI and +PI (radians)
                // https://developer.apple.com/documentation/accelerate/1450132-vdsp_zvphasd/
                vDSP_zvphasD(&tempSplitComplex, 1, &phase, 1, halfN_);
            }
        }
        
        // ----------------------------------------------------------------
        // Bandpass Filtering
        // ----------------------------------------------------------------
        
        // Get the Frequencies for the current samplerate
        let freqs = getFrequencies(N, samplerate: samplerate)
        // Get a Bandpass Filter
        let bandPassFilter = generateBandPassFilter(freqs)
        
        // Multiply phase and magnitude by the bandpass filter
        mag = mul(mag, y: bandPassFilter.0)
        phase = mul(phase, y: bandPassFilter.0)
        
        // Output Variables
        let filteredSpectrum = mul(fullSpectrum, y: bandPassFilter.0)
        let filteredPhase = phase
        
        // ----------------------------------------------------------------
        // Determine Maximum Frequency
        // ----------------------------------------------------------------
        let (maxAmplitude, maxIndex) = max(filteredSpectrum)
        let maxFrequency = freqs[maxIndex]
        let maxPhase = filteredPhase[maxIndex]
        
        let phasePercent = abs(((maxPhase / .pi) + 0.5) * 100.0) // Normalize to 0-100%
        
        print("Amplitude: \(maxAmplitude)")
        print("Frequency: \(maxFrequency)")
        print("Phase: \(String(format: "%.2f", phasePercent))%")
        
        vDSP_DFT_DestroySetupD(dftForwardSetup)
        
        
        // ----------------------------------------------------------------
        // Convert from polar coordinates back to rectangular coordinates.
        // ----------------------------------------------------------------
        
        mag.withUnsafeMutableBufferPointer { magPtr in
            phase.withUnsafeMutableBufferPointer { phasePtr in
                var tempSplitComplex = SplitComplex(realp: magPtr.baseAddress!,
                                                    imagp: phasePtr.baseAddress!)
                
                var complexAsValue : UnsafeMutablePointer<Sample>? = nil
                
                tempComplex.withUnsafeMutableBytes {
                    complexAsValue = $0.baseAddress?.bindMemory(to: Sample.self, capacity: N)
                }
                
                vDSP_ztocD(&tempSplitComplex, 1, &tempComplex, 2, halfN_);
                vDSP_rectD(complexAsValue!, 2, complexAsValue!, 2, halfN_);
                vDSP_ctozD(&tempComplex, 2, &tempSplitComplex, 1, halfN_);
            }
        }
        
        // ----------------------------------------------------------------
        // Do Inverse DFT
        // ----------------------------------------------------------------
        
        // Do complex->real inverse out-of-place DFT.
        vDSP_DFT_ExecuteD(dftInverseSetup, mag, phase, &tempReal, &tempImag)
        
        // Create result
        var result: [Sample] = [Sample](repeating: 0.0, count: N)
       
        // The Accelerate DFT leaves the result’s even and odd values split.
        // Here we zip them together into a real vector.
        // https://developer.apple.com/documentation/accelerate/performing_fourier_transforms_on_interleaved-complex_data
        tempReal.withUnsafeMutableBufferPointer { tempRealPtr in
            tempImag.withUnsafeMutableBufferPointer { tempImagPtr in
                
                var tempSplitComplex =
                    SplitComplex(realp: tempRealPtr.baseAddress!,
                                 imagp: tempImagPtr.baseAddress!)
                
                var resultAsComplex: UnsafeMutablePointer<Complex>? = nil
                
                result.withUnsafeMutableBytes {
                    resultAsComplex =
                        $0.baseAddress?.bindMemory(to: Complex.self,
                                                   capacity: N)
                }
                
                vDSP_ztocD(&tempSplitComplex, vDSP_Stride(1),
                           resultAsComplex!, vDSP_Stride(2),
                           halfN_);
            }
        }
        
        // Neither the forward nor inverse DFT does any scaling. We compensate for that here.
        var scale = 0.5/Sample(N);
        var resultScaled = [Sample](repeating: 0.0, count: N);
        vDSP_vsmulD(&result, 1, &scale, &resultScaled, 1, N_);
        result = resultScaled
        
        #if false
        // Print Result
        for k in 0 ..< N {
            print("\(k)   \(inRealValues[k])     \(result[k])")
        }
        #endif
        
        vDSP_DFT_DestroySetupD(dftInverseSetup)
    }
        
    
    // The bandpass frequencies
    let lowerFreq : Sample = 3.0
    let higherFreq: Sample = 5.0
    
    // Some Math functions on Arrays
    func mul(_ x: [Sample], y: [Sample]) -> [Sample] {
        var results = [Sample](repeating: 0.0, count: x.count)
        vDSP_vmulD(x, 1, y, 1, &results, 1, vDSP_Length(x.count))
        
        return results
    }
    
    func sqrt(_ x: [Sample]) -> [Sample] {
        var results = [Sample](repeating: 0.0, count: x.count)
        vvsqrt(&results, x, [Int32(x.count)])
        
        return results
    }
    
    func max(_ x: [Sample]) -> (Sample, Int) {
        var result: Sample = 0.0
        var idx : vDSP_Length = vDSP_Length(0)
        vDSP_maxviD(x, 1, &result, &idx, vDSP_Length(x.count))
        
        return (result, Int(idx))
    }
}
