using System;

namespace AASharp
{
    /// <summary>
    /// This class provides for the calculation of sidereal time. This refers to Chapter 12 in the book.
    /// </summary>
    public static class AASSidereal
    {
        /// <param name="JD">The Julian Day in Universal time to calculate for.</param>
        /// <returns>The Mean Greenwich Sidereal Time, that is, the Greenwich hour angle of the mean vernal point (the intersection of the ecliptic of the date with the mean equator of the date), expressed in hours.</returns>
        public static double MeanGreenwichSiderealTime(double JD)
        {
            //Get the Julian day for the same day at midnight
            long Year = 0;
            long Month = 0;
            long Day = 0;
            long Hour = 0;
            long Minute = 0;
            double Second = 0;

            AASDate date = new AASDate();
            date.Set(JD, AASDate.AfterPapalReform(JD));
            date.Get(ref Year, ref Month, ref Day, ref Hour, ref Minute, ref Second);
            date.Set(Year, Month, Day, 0, 0, 0, date.InGregorianCalendar);
            double JDMidnight = date.Julian;

            //Calculate the sidereal time at midnight
            double T = (JDMidnight - 2451545) / 36525;
            double TSquared = T * T;
            double TCubed = TSquared * T;
            double Value = 100.46061837 + (36000.770053608 * T) + (0.000387933 * TSquared) - (TCubed / 38710000);

            //Adjust by the time of day
            Value += (((Hour * 15) + (Minute * 0.25) + (Second * 0.0041666666666666666666666666666667)) * 1.00273790935);

            Value = AASCoordinateTransformation.DegreesToHours(Value);

            return AASCoordinateTransformation.MapTo0To24Range(Value);
        }

        /// <param name="JD">The Julian Day in Universal time to calculate for.</param>
        /// <returns>The Apparent Greenwich Sidereal Time, that is, the Greenwich hour angle of the true vernal equinox, expressed in hours.</returns>
        public static double ApparentGreenwichSiderealTime(double JD)
        {
            double MeanObliquity = AASNutation.MeanObliquityOfEcliptic(JD);
            double TrueObliquity = MeanObliquity + AASNutation.NutationInObliquity(JD) / 3600;
            double NutationInLongitude = AASNutation.NutationInLongitude(JD);

            double Value = MeanGreenwichSiderealTime(JD) + (NutationInLongitude * Math.Cos(AASCoordinateTransformation.DegreesToRadians(TrueObliquity)) / 54000);
            return AASCoordinateTransformation.MapTo0To24Range(Value);
        }
    }
}
