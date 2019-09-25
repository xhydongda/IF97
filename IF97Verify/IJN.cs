using System;
using System.Collections.Generic;
using System.Text;

namespace IF97
{
    /// <summary>
    /// Structure for the double indexed state equation coefficients.
    /// </summary>
    public struct IJn
    {
        public IJn(int i, int j, double n)
        {
            I = i;
            J = j;
            this.n = n;
        }
        /// <summary>
        /// The first index.
        /// </summary>
        public int I;
        /// <summary>
        /// The second index.
        /// </summary>
        public int J;
        /// <summary>
        /// The leading numerical constant.
        /// </summary>
        public double n;
    }
}
