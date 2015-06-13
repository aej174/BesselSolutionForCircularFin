//
//  ViewController.swift
//  BesselSolutionForCircularFin
//
//  Created by Allan Jones on 6/12/15.
//  Copyright (c) 2015 Allan Jones. All rights reserved.
//

import UIKit

class ViewController: UIViewController {
  
  @IBOutlet weak var finODTextField: UITextField!
  @IBOutlet weak var finIDTextField: UITextField!
  @IBOutlet weak var finThicknessTextField: UITextField!
  @IBOutlet weak var thermalConductivityTextField: UITextField!
  @IBOutlet weak var heatTransferCoefficientTextField: UITextField!
  
  @IBOutlet weak var hotFluidTemperatureTextField: UITextField!
  @IBOutlet weak var ambientTemperatureTextField: UITextField!
  
  @IBOutlet weak var finTipTemperatureTextField: UITextField!
  @IBOutlet weak var heatTransferRateTextField: UITextField!
  
  let pi = 3.1415926
  
  let numberOfCalcs = 10
  
  var i0: Double = 0.0
  var i1: Double = 0.0
  var k0: Double = 0.0
  var k1: Double = 0.0
  
  override func viewDidLoad() {
    super.viewDidLoad()
    // Do any additional setup after loading the view, typically from a nib.
  }

  override func didReceiveMemoryWarning() {
    super.didReceiveMemoryWarning()
    // Dispose of any resources that can be recreated.
  }
  
  @IBAction func calculateButtonPressed(sender: UIButton) {
    
    let finOD = Double((finODTextField.text as NSString).doubleValue)
    let finID = Double((finIDTextField.text as NSString).doubleValue)
    let thicknessOfFin = Double((finThicknessTextField.text as NSString).doubleValue)
    let thermalConductivity = Double((thermalConductivityTextField.text as NSString).doubleValue)
    let heatTransferCoefficient = Double((heatTransferCoefficientTextField.text as NSString).doubleValue)
    let hotFluidTemp = Double((hotFluidTemperatureTextField.text as NSString).doubleValue)
    let ambientTemp = Double((ambientTemperatureTextField.text as NSString).doubleValue)
    
    let r0 = finID / (2.0 * 12.0)
    let r1 = finOD / (2.0 * 12.0)
    let finThickness = thicknessOfFin / 12.0
    
    let beta2 = 2.0 * heatTransferCoefficient / (finThickness * thermalConductivity)
    let beta = sqrt(beta2)       //dimension:  ft^-1
    println("r0 = \(r0), r1 = \(r1), beta = \(beta)")
    
    let x0 = beta * r0
    let x1 = beta * r1
    
    var i0r0 = computeI0(x0)
    var i1r0 = computeI1(x0)
    var k0r0 = computeK0(x0)
    var k1r0 = computeK1(x0)
    
    var i0r1 = computeI0(x1)
    var i1r1 = computeI1(x1)
    var k0r1 = computeK0(x1)
    var k1r1 = computeK1(x1)
    
    println("x0 = \(x0), x1 = \(x1)")
    println("i0r0 = \(i0r0), i1r0 = \(i1r0), k0r0 = \(k0r0), k1r0 = \(k1r0)")
    println("i0r1 = \(i0r1), i1r1 = \(i1r1), k0r1 = \(k0r1), k1r1 = \(k1r1)")
    
    let denom = (k1r1 * i0r0) + (i1r1 * k0r0)
    let c1 = k1r1 / denom
    let c2 = i1r1 / denom
    
    println("denom = \(denom), c1 = \(c1), c2 = \(c2)")
    
    let deltaR = (r1 - r0) / Double(numberOfCalcs)
    
    for var i = 0; i < numberOfCalcs; ++i {
      var r = r0 + Double(i) * deltaR
      var rx = beta * r
      var i0rx = computeI0(rx)
      var k0rx = computeK0(rx)
      var theta = c1 * i0rx + c2 * k0rx
      var temp = ambientTemp + theta * (hotFluidTemp - ambientTemp)
      println("r = \(r), temp = \(temp)")
    }
    var finTipTempRatio = c1 * i0r1 + c2 * k0r1
    var finTipTemp = ambientTemp + finTipTempRatio * (hotFluidTemp - ambientTemp)
    finTipTemperatureTextField.text = "\(finTipTemp)"
    println("r1 = \(r1), temp = \(finTipTemp)")
    
    var heatFluxAtBase = thermalConductivity * (hotFluidTemp - ambientTemp) * beta * (i1r1 * k1r0 - k1r1 * i1r0) / denom
    let finBaseArea = 2.0 * pi * r0 * finThickness
    println("heat flux at base = \(heatFluxAtBase), BTU/hr-ft2")
    println("fin base area = \(finBaseArea) ft2")
    let heatTransferRate = finBaseArea * heatFluxAtBase
    heatTransferRateTextField.text = "\(heatTransferRate)"
  }
  
  //MARK: Bessel Functions
  
  func computeI0(xx: Double) -> Double {
    let tt = xx / 3.75
    
    if (xx > -3.75) && (xx < 3.75) {
      let c1 = 3.5156229
      let c2 = 3.0899424
      let c3 = 1.2067492
      let c4 = 0.2659732
      let c5 = 0.0360768
      let c6 = 0.0045813
      
      i0 = 1.0 + c1 * pow(tt, 2.0) + c2 * pow(tt, 4.0) + c3 * pow(tt, 6.0) + c4 * pow(tt, 8.0) + c5 * pow(tt, 10.0) + c6 * pow(tt, 12.0)
    }
    else {
      let c1 = 0.39894228
      let c2 = 0.01328592
      let c3 = 0.00225319
      let c4 = -0.00157565
      let c5 = 0.00916281
      let c6 = -0.02057706
      let c7 = 0.02635537
      let c8 = -0.01647633
      let c9 = 0.00392377
      
      let c10 = sqrt(xx) * exp(-xx)
      let sum = c1 + c2 / tt + c3 / pow(tt, 2.0) + c4 / pow(tt, 3.0) + c5 / pow(tt, 4.0) + c6 / pow(tt, 5.0) + c7 / pow(tt, 6.0) + c8 / pow(tt, 7.0) + c9 / pow(tt, 8.0)
      i0 = sum / c10
    }
    return i0
  }
  
  func computeI1(xx: Double) -> Double {
    let tt = xx / 3.75
    
    if (xx > -3.75) && (xx < 3.75) {
      let d1 = 0.87890594
      let d2 = 0.51498869
      let d3 = 0.15084934
      let d4 = 0.02658733
      let d5 = 0.00301532
      let d6 = 0.00032411
      
      let sum = 0.5 + d1 * pow(tt, 2.0) + d2 * pow(tt, 4.0) + d3 * pow(tt, 6.0) + d4 * pow(tt, 8.0) + d5 * pow(tt, 10.0) + d6 * pow(tt, 12.0)
      i1 = sum * xx
    }
    else {
      let d1 = 0.39894228
      let d2 = -0.03988024
      let d3 = -0.00362018
      let d4 = 0.00163801
      let d5 = -0.01031555
      let d6 = 0.02282967
      let d7 = -0.02895312
      let d8 = 0.01787654
      let d9 = -0.00420059
      
      let d10 = sqrt(xx) * exp(-xx)
      let sum = d1 + d2 / tt + d3 / pow(tt, 2.0) + d4 / pow(tt, 3.0) + d5 / pow(tt, 4.0) + d6 / pow(tt, 5.0) + d7 / pow(tt, 6.0) + d8 / pow(tt, 7.0) + d9 / pow(tt, 8.0)
      i1 = sum / d10
    }
    return i1
  }
  
  func computeK0(xx: Double) -> Double {
    let pp = xx / 2.0
    
    if xx < 2.0 {
      let f1 = -0.57721566
      let f2 = 0.42278420
      let f3 = 0.23069756
      let f4 = 0.03488590
      let f5 = 0.00262698
      let f6 = 0.00010750
      let f7 = 0.00000740
      
      let f8 = -log(pp) * i0
      k0 = f8 + f1 + f2 * pow(pp, 2.0) + f3 * pow(pp, 4.0) + f4 * pow(pp, 6.0) + f5 * pow(pp, 8.0) + f6 * pow(pp, 10.0) + f7 * pow(pp, 12.0)
    }
    else {
      let f1 = 1.25331414
      let f2 = -0.07832358
      let f3 = 0.02189568
      let f4 = -0.01062446
      let f5 = 0.00587872
      let f6 = -0.00251540
      let f7 = 0.00053208
      
      let f8 = sqrt(xx) * exp(xx)
      let sum = f1 + f2 / pp + f3 / pow(pp, 2.0) + f4 / pow(pp, 3.0) + f5 / pow(pp, 4.0) + f6 / pow(pp, 5.0) + f7 / pow(pp, 6.0)
      k0 = sum / f8
    }
    return k0
  }
  
  func computeK1(xx: Double) -> Double {
    let pp = xx / 2.0
    
    if (xx < 2.0) {
      let g1 = 0.15443144
      let g2 = -0.67278579
      let g3 = -0.18156897
      let g4 = -0.01919402
      let g5 = -0.00110404
      let g6 = -0.00004686
      
      let g7 = xx * log(pp) * i1
      let sum = g7 + 1.0 + g1 * pow(pp, 2.0) + g2 * pow(pp, 4.0) + g3 * pow(pp, 6.0) + g4 * pow(pp, 8.0) + g5 * pow(pp, 10.0) + g6 * pow(pp, 12.0)
      k1 = sum / xx
    }
    else {
      let g1 = 1.25331414
      let g2 = 0.23498619
      let g3 = -0.03655620
      let g4 = 0.01504268
      let g5 = -0.00780353
      let g6 = 0.00325614
      let g7 = -0.00068245
      
      let g8 = sqrt(xx) * exp(xx)
      let sum = g1 + g2 / pp + g3 / pow(pp, 2.0) + g4 / pow(pp, 3.0) + g5 / pow(pp, 4.0) + g6 / pow(pp, 5.0) + g7 / pow(pp, 6.0)
      k1 = sum / g8
    }
    return k1
  }

}

