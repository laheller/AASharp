using System;

namespace AASharp {
    /// <summary>
    /// This class provides for calculation of the position of an object in an orbit which is hyperbolic i.e. with an eccentricity greater than one. One example being the 2025 interstellar comet 3I/ATLAS.
    /// </summary>
    public static class AASHyperbolic {
        private static double CalculateKeplers(double M, double e, double epsilon = 0.000001) {
            var bRecalc = true;
            var Hk = M;
            while (bRecalc) {
                var Hk1 = Hk - ((e * Math.Sinh(Hk) - Hk - M) / ((e * Math.Cosh(Hk)) - 1));

                //Prepare for the next loop around
                bRecalc = (Math.Abs(Hk1 - Hk) > epsilon);
                Hk = Hk1;
            }

            return Hk;
        }

        private static void CalculateTrueAnomalyAndRadius(double JD, ref AASHyperbolicObjectElements elements, ref double v, ref double r, double epsilon = 0.000001) {
            var a = elements.q / (elements.e - 1);
            var n = Math.Sqrt(0.0002959122082855911025 / Math.Pow(Math.Abs(a), 3));
            var M = (JD - elements.T) * n;
            var H = CalculateKeplers(M, elements.e, epsilon);
            v = 2 * Math.Atan(Math.Sqrt((elements.e + 1) / (elements.e - 1)) * Math.Tanh(H / 2));
            r = Math.Abs(a) * (elements.e * Math.Cosh(H) - 1);
        }

        /// <param name="JD">The date in Dynamical time to calculate for.</param>
        /// <param name="elements">An instance of the AASHyperbolicObjectElements class with the orbit elements.</param>
        /// <param name="bHighPrecision">If true then use the full VSOP87 theory instead of the truncated version as provided in Meeus's book.</param>
        /// <param name="epsilon">The epsilon value passed to the CalculateKeplers method.</param>
        /// <returns></returns>
        public static AASHyperbolicObjectDetails Calculate(double JD, ref AASHyperbolicObjectElements elements, bool bHighPrecision, double epsilon = 0.000001) {
            var Epsilon = AASNutation.MeanObliquityOfEcliptic(elements.JDEquinox);
            var JD0 = JD;

            //What will be the return value
            var details = new AASHyperbolicObjectDetails();

            Epsilon = AASCoordinateTransformation.DegreesToRadians(Epsilon);
            var omega = AASCoordinateTransformation.DegreesToRadians(elements.omega);
            var w = AASCoordinateTransformation.DegreesToRadians(elements.w);
            var i = AASCoordinateTransformation.DegreesToRadians(elements.i);

            var sinEpsilon = Math.Sin(Epsilon);
            var cosEpsilon = Math.Cos(Epsilon);
            var sinOmega = Math.Sin(omega);
            var cosOmega = Math.Cos(omega);
            var cosi = Math.Cos(i);
            var sini = Math.Sin(i);

            var F = cosOmega;
            var G =  sinOmega * cosEpsilon;
            var H = sinOmega * sinEpsilon;
            var P = -sinOmega * cosi;
            var Q = (cosOmega * cosi * cosEpsilon) - (sini * sinEpsilon);
            var R = (cosOmega * cosi * sinEpsilon) + (sini * cosEpsilon);
            var a = Math.Sqrt((F * F) + (P * P));
            var b = Math.Sqrt((G * G) + (Q * Q));
            var c = Math.Sqrt((H * H) + (R * R));
            var A = Math.Atan2(F, P);
            var B = Math.Atan2(G, Q);
            var C = Math.Atan2(H, R);

            var SunCoord = AASSun.EquatorialRectangularCoordinatesAnyEquinox(JD, elements.JDEquinox, bHighPrecision);
            for (int j = 0; j < 2; j++) {
                var v = 0.0;
                var r = 0.0;
                CalculateTrueAnomalyAndRadius(JD0, ref elements, ref v, ref r, epsilon);
                var x = r * a * Math.Sin(A + w + v);
                var y = r * b * Math.Sin(B + w + v);
                var z = r * c * Math.Sin(C + w + v);

                if (j == 0) {
                    details.HeliocentricRectangularEquatorial.X = x;
                    details.HeliocentricRectangularEquatorial.Y = y;
                    details.HeliocentricRectangularEquatorial.Z = z;

                    //Calculate the heliocentric ecliptic coordinates also
                    var u = w + v;
                    var cosu = Math.Cos(u);
                    var sinu = Math.Sin(u);

                    details.HeliocentricRectangularEcliptical.X = r * ((cosOmega * cosu) - (sinOmega * sinu * cosi));
                    details.HeliocentricRectangularEcliptical.Y = r * ((sinOmega * cosu) + (cosOmega * sinu * cosi));
                    details.HeliocentricRectangularEcliptical.Z = r * sini * sinu;

                    details.HeliocentricEclipticLongitude = AASCoordinateTransformation.MapTo0To360Range(AASCoordinateTransformation.RadiansToDegrees(Math.Atan2(details.HeliocentricRectangularEcliptical.Y, details.HeliocentricRectangularEcliptical.X)));
                    details.HeliocentricEclipticLatitude = AASCoordinateTransformation.RadiansToDegrees(Math.Asin(details.HeliocentricRectangularEcliptical.Z / r));
                }

                var psi = SunCoord.X + x;
                var psi2 = psi * psi;
                var nu = SunCoord.Y + y;
                var nu2 = nu * nu;
                var sigma = SunCoord.Z + z;

                var Alpha = Math.Atan2(nu, psi);
                Alpha = AASCoordinateTransformation.RadiansToDegrees(Alpha);
                var Delta = Math.Atan2(sigma, Math.Sqrt(psi2 + nu2));
                Delta = AASCoordinateTransformation.RadiansToDegrees(Delta);
                var Distance = Math.Sqrt(psi2 + nu2 + (sigma * sigma));
                var Distance2 = Distance * Distance;

                if (j == 0) {
                    details.TrueGeocentricRA = AASCoordinateTransformation.MapTo0To24Range(Alpha / 15);
                    details.TrueGeocentricDeclination = Delta;
                    details.TrueGeocentricDistance = Distance;
                    details.TrueGeocentricLightTime = AASElliptical.DistanceToLightTime(Distance);
                }
                else {
                    details.AstrometricGeocentricRA = AASCoordinateTransformation.MapTo0To24Range(Alpha / 15);
                    details.AstrometricGeocentricDeclination = Delta;
                    details.AstrometricGeocentricDistance = Distance;
                    details.AstrometricGeocentricLightTime = AASElliptical.DistanceToLightTime(Distance);

                    var RES = Math.Sqrt((SunCoord.X * SunCoord.X) + (SunCoord.Y * SunCoord.Y) + (SunCoord.Z * SunCoord.Z));
                    var RES2 = RES * RES;
                    var r2 = r * r;

                    details.Elongation = AASCoordinateTransformation.RadiansToDegrees(Math.Acos((RES2 + Distance2 - r2) / (2 * RES * Distance)));
                    details.PhaseAngle = AASCoordinateTransformation.RadiansToDegrees(Math.Acos((r2 + Distance2 - RES2) / (2 * r * Distance)));
                }

                if (j == 0) //Prepare for the next loop around
                    JD0 = JD - details.TrueGeocentricLightTime;
            }

            return details;
        }
    }
}
