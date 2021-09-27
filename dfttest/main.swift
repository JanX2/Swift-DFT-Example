//
//  main.swift
//  dfttest
//
//  Created by Christopher Helf on 17.08.15.
//  Copyright (c) 2015-2021 Christopher Helf. All rights reserved.
//  Copyright (c) 2021-Present Jan Weiß. All rights reserved.
//

import Foundation
import Accelerate

var dft = DFT()

let n = 512 // Should be a power of two for best performance.
let frequency1 = 4.0
let phase1 = 0.0
let amplitude1 = 8.0
let seconds = 2.0
let rate = Double(n)/seconds

var sineWave = (0..<n).map {
    amplitude1 * sin(2.0 * .pi / rate * DFT.Sample($0) * frequency1 + phase1)
}

dft.calculate(sineWave, samplerate: rate)

