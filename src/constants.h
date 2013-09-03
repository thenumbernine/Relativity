#pragma once

//universal constants

const double speedOfLightInMPerS = 299792458.;
const double gravitationalConstantInM3PerKgS2 = 6.67384e-11;

//conversions to meters
const double metersPerS = speedOfLightInMPerS;
const double metersPerKg = gravitationalConstantInM3PerKgS2 / (metersPerS * metersPerS);

//properties of the sun
const double sunMassInKg = 1.989e+30;
const double sunMassInM = sunMassInKg * metersPerKg;
const double sunRadiusInM = 6.955e+8;
const double sunVolumeInM3 = 4. / 3. * M_PI * sunRadiusInM * sunRadiusInM * sunRadiusInM;
const double sunDensityInM_2 = sunMassInKg * metersPerKg / sunVolumeInM3;

const double sunRotationInRad_S = 2.8e-6;
const double sunRotationInM = sunRotationInRad_S / metersPerS;
const double sunAngularMomentumInKgM2_S = .4 * sunMassInKg * sunRadiusInM * sunRadiusInM * sunRotationInRad_S;
const double sunAngularMomentumInM = .4 * sunMassInM * sunRadiusInM * sunRadiusInM * sunRotationInM;

