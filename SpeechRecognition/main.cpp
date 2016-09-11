//
//  main.cpp
//  SpeechRecognition
//
//  Created by hty on 9/2/16.
//  Copyright © 2016 hty. All rights reserved.
//

#include <iostream>
#include "portaudio.h"
#include "captureData.h"
#include "fft.h"
#include "featureExtraction.h"

using namespace std;


int main(int argc, const char * argv[]) {
//    capture();
//    test();
    featureExtraction();
}
