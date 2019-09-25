namespace IF97
{
    public enum IF97Parameters
    {
        /// <summary>
        /// 密度ρ.
        /// </summary>
        d,
        /// <summary>
        /// 比焓.
        /// </summary>
        h,
        /// <summary>
        /// 绝对温度.
        /// </summary>
        T,
        /// <summary>
        /// 压力.
        /// </summary>
        p,
        /// <summary>
        /// 比熵.
        /// </summary>
        s,
        /// <summary>
        /// 比内能.
        /// </summary>
        u,
        /// <summary>
        /// 比体积.
        /// </summary>
        v,
        /// <summary>
        /// 等压比热.
        /// </summary>
        cp,
        /// <summary>
        /// 等容比热.
        /// </summary>
        cv,
        /// <summary>
        /// 声速.
        /// </summary>
        w,
        drhodp,
        /// <summary>
        /// 动力粘度,μ.
        /// </summary>
        mu,
        /// <summary>
        /// 热导率.
        /// </summary>
        tc
    }

    /// <summary>
    /// Saturated Liquid/Vapor state determination
    /// </summary>
    public enum SatState
    {
        NONE, LIQUID, VAPOR
    }
}
