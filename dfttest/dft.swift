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
    
    fileprivate func getFrequencies(_ N: Int, fps: Double) -> [Double] {
        // Create an Array with the Frequencies
        let freqs = (0..<N/2).map {
            fps/Double(N)*Double($0)
        }
        
        return freqs
    }
    
    fileprivate func generateBandPassFilter(_ freqs: [Double]) -> ([Double], Int, Int) {
        var minIdx = freqs.count+1
        var maxIdx = -1
        
        let bandPassFilter: [Double] = freqs.map {
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
    
    func calculate(_ inRealValues: [Double], fps: Double) {
        // ----------------------------------------------------------------
        // Size Variables
        // ----------------------------------------------------------------
        let N = inRealValues.count
        let halfN = N/2
        let N_ = vDSP_Length(N)
        let halfN_ = vDSP_Length(halfN)
        
        // --------------------------------------------------------------------
        // Split even and odd indexes of our input values into their own arrays
        // --------------------------------------------------------------------
        var inEven = [Double](repeating: 0.0, count:halfN)
        var inOdd = [Double](repeating: 0.0, count:halfN)
        inEven.withUnsafeMutableBufferPointer { inEvenPtr in
            inOdd.withUnsafeMutableBufferPointer { inOddPtr in
                
                var tempSplitComplex =
                    DSPDoubleSplitComplex(realp: inEvenPtr.baseAddress!,
                                          imagp: inOddPtr.baseAddress!)
                
                var valuesAsComplex: UnsafePointer<DSPDoubleComplex>? = nil
                
                inRealValues.withUnsafeBytes {
                    valuesAsComplex =
                        $0.baseAddress?.bindMemory(to: DSPDoubleComplex.self,
                                                   capacity: N)
                }
                
                vDSP_ctozD(valuesAsComplex!, vDSP_Stride(2),
                           &tempSplitComplex, vDSP_Stride(1),
                           halfN_)
                
            }
        }
        
        // ----------------------------------------------------------------
        // DFT & Variables Setup
        // ----------------------------------------------------------------
        let dftForwardSetup: vDSP_DFT_Setup = vDSP_DFT_zrop_CreateSetupD(nil, N_, .FORWARD)!
        let dftInverseSetup: vDSP_DFT_Setup = vDSP_DFT_zrop_CreateSetupD(dftForwardSetup, N_, .INVERSE)!
        
        // We need the complex buffer in two different data layouts!
        var tempComplex : [DSPDoubleComplex] = [DSPDoubleComplex](repeating: DSPDoubleComplex(), count: halfN)
        
        var tempReal : [Double] = [Double](repeating: 0.0, count: halfN)
        var tempImag : [Double] = [Double](repeating: 0.0, count: halfN)
        
        // For polar coordinates
        var mag : [Double] = [Double](repeating: 0.0, count: halfN)
        var phase : [Double] = [Double](repeating: 0.0, count: halfN)
        
        // ----------------------------------------------------------------
        // Forward DFT
        // ----------------------------------------------------------------
        
        // Do real->complex forward DFT
        vDSP_DFT_ExecuteD(dftForwardSetup, &inEven, &inOdd, &tempReal, &tempImag)
        
        // ----------------------------------------------------------------
        // Get the Frequency Spectrum
        // ----------------------------------------------------------------
        var fullSpectrum = [Double](repeating: 0.0, count: halfN)
        
        tempReal.withUnsafeMutableBufferPointer { tempRealPtr in
            tempImag.withUnsafeMutableBufferPointer { tempImagPtr in
                
                var tempSplitComplex = DSPDoubleSplitComplex(realp: tempRealPtr.baseAddress!,
                                                             imagp: tempImagPtr.baseAddress!)
                
                var dftMagnitudes = [Double](repeating: 0.0, count: halfN)
                vDSP_zvmagsD(&tempSplitComplex, 1, &dftMagnitudes, 1, halfN_);
                
                // vDSP_zvmagsD returns squares of the DFT magnitudes, so take the root here
                let roots = sqrt(dftMagnitudes)
                
                // Normalize the Amplitudes
                vDSP_vsmulD(roots, vDSP_Stride(1), [1.0 / Double(N)], &fullSpectrum, 1, halfN_)
                
                // ----------------------------------------------------------------
                // Convert from complex/rectangular (real, imaginary) coordinates
                // to polar (magnitude and phase) coordinates.
                // ----------------------------------------------------------------
                
                vDSP_zvabsD(&tempSplitComplex, vDSP_Stride(1),
                            &mag, vDSP_Stride(1),
                            halfN_);
                
                // Beware: Outputted phase here between -PI and +PI
                // https://developer.apple.com/library/prerelease/ios/documentation/Accelerate/Reference/vDSPRef/index.html#//apple_ref/c/func/vDSP_zvphasD
                vDSP_zvphasD(&tempSplitComplex, 1, &phase, 1, halfN_);
            }
        }
        
        // ----------------------------------------------------------------
        // Bandpass Filtering
        // ----------------------------------------------------------------
        
        // Get the Frequencies for the current Framerate
        let freqs = getFrequencies(N, fps: fps)
        // Get a Bandpass Filter
        let bandPassFilter = generateBandPassFilter(freqs)
        
        // Multiply phase and magnitude with the bandpass filter
        mag = mul(mag, y: bandPassFilter.0)
        phase = mul(phase, y: bandPassFilter.0)
        
        // Output Variables
        let filteredSpectrum = mul(fullSpectrum, y: bandPassFilter.0)
        let filteredPhase = phase
        
        // ----------------------------------------------------------------
        // Determine Maximum Frequency
        // ----------------------------------------------------------------
        let maxFrequencyResult = max(filteredSpectrum)
        let maxFrequency = freqs[maxFrequencyResult.1]
        let maxPhase = filteredPhase[maxFrequencyResult.1]
        
        print("Amplitude: \(maxFrequencyResult.0)")
        print("Frequency: \(maxFrequency)")
        print("Phase: \(maxPhase + .pi / 2)")
        
        vDSP_DFT_DestroySetupD(dftForwardSetup)
        
        
        mag.withUnsafeMutableBufferPointer { magPtr in
            phase.withUnsafeMutableBufferPointer { phasePtr in
                // ----------------------------------------------------------------
                // Convert from polar coordinates back to rectangular coordinates.
                // ----------------------------------------------------------------
                
                var tempSplitComplex = DSPDoubleSplitComplex(realp: magPtr.baseAddress!,
                                                             imagp: phasePtr.baseAddress!)
                
                var complexAsValue : UnsafeMutablePointer<Double>? = nil
                
                tempComplex.withUnsafeMutableBytes {
                    complexAsValue = $0.baseAddress?.bindMemory(to: Double.self, capacity: N)
                }
                
                vDSP_ztocD(&tempSplitComplex, 1, &tempComplex, 2, halfN_);
                vDSP_rectD(complexAsValue!, 2, complexAsValue!, 2, halfN_);
                vDSP_ctozD(&tempComplex, 2, &tempSplitComplex, 1, halfN_);
            }
        }
        
        // ----------------------------------------------------------------
        // Do Inverse DFT
        // ----------------------------------------------------------------
        var outEven = [Double](repeating: 0.0, count:halfN) // Re-use tempReal and tempImag?
        var outOdd = [Double](repeating: 0.0, count:halfN)
        
        // Do complex->real inverse DFT.
        vDSP_DFT_ExecuteD(dftInverseSetup, &mag, &phase, &outEven, &outOdd)
        
        // Create result
        var result: [Double] = [Double](repeating: 0.0, count: N)
       
        // The Accelerate DFT leaves the result’s even and odd values split. Here we zip them together into a real vector.
        outEven.withUnsafeMutableBufferPointer { outEvenPtr in
            outOdd.withUnsafeMutableBufferPointer { outOddPtr in
                
                var tempSplitComplex =
                    DSPDoubleSplitComplex(realp: outEvenPtr.baseAddress!,
                                          imagp: outOddPtr.baseAddress!)
                
                var resultAsComplex: UnsafeMutablePointer<DSPDoubleComplex>? = nil
                
                result.withUnsafeMutableBytes {
                    resultAsComplex =
                        $0.baseAddress?.bindMemory(to: DSPDoubleComplex.self,
                                                                 capacity: N)
                }
                
                vDSP_ztocD(&tempSplitComplex, vDSP_Stride(1),
                           resultAsComplex!, vDSP_Stride(2),
                           halfN_);
            }
        }
        
        // Neither the forward nor inverse DFT does any scaling. Here we compensate for that.
        var scale : Double = 0.5/Double(N);
        var resultScaled = [Double](repeating: 0.0, count: N);
        vDSP_vsmulD(&result, 1, &scale, &resultScaled, 1, N_);
        result = resultScaled
        
        // Print Result
        for k in 0 ..< N {
            print("\(k)   \(inRealValues[k])     \(result[k])")
        }
        
        vDSP_DFT_DestroySetupD(dftInverseSetup)
    }
        
    
    // The bandpass frequencies
    let lowerFreq : Double = 3
    let higherFreq: Double = 5
    
    // Some Math functions on Arrays
    func mul(_ x: [Double], y: [Double]) -> [Double] {
        var results = [Double](repeating: 0.0, count: x.count)
        vDSP_vmulD(x, 1, y, 1, &results, 1, vDSP_Length(x.count))
        
        return results
    }
    
    func sqrt(_ x: [Double]) -> [Double] {
        var results = [Double](repeating: 0.0, count: x.count)
        vvsqrt(&results, x, [Int32(x.count)])
        
        return results
    }
    
    func max(_ x: [Double]) -> (Double, Int) {
        var result: Double = 0.0
        var idx : vDSP_Length = vDSP_Length(0)
        vDSP_maxviD(x, 1, &result, &idx, vDSP_Length(x.count))
        
        return (result, Int(idx))
    }
}
