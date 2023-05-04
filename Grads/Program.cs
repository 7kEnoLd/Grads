using System.Runtime.CompilerServices;
using System.Xml.XPath;

public class Program
{
    #region 底层方法
    static int scope = 100;

    static List<double> fig = new List<double>();

    static int iter = 100000;

    static double tol = 0.001;

    public static double GetFuncNum(int scope, double[] x)
    {
        double f = 0;
        for (int i = 0; i < scope; i++)
        {
            f += Math.Pow(-x[i] + x[i + 1] + x[i + 2], 2) + Math.Pow(x[i] - x[i + 1] + x[i + 2], 2) + Math.Pow(x[i] + x[i + 1] - x[i + 2], 2);
        }
        return f;
    }

    public static double[] GetFuncDen(int scope, double[] x)
    {
        double[] f = new double[scope + 2];
        for (int i = 0; i < scope; i++)
        {
            f[i] += 6 * x[i] - 2 * x[i + 1] - 2 * x[i + 2];
            f[i + 1] += -2 * x[i] + 6 * x[i + 1] - 2 * x[i + 2];
            f[i + 2] += -2 * x[i] - 2 * x[i + 1] + 6 * x[i + 2];
        }
        return f;
    }

    static double[] GetX(int scope)
    {
        double[] x = new double[scope + 2];
        if (scope % 3 == 1)
        {
            int count = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                switch (count)
                {
                    case 0:
                        x[i] = 1; break;
                    case 1:
                        x[i] = 2; break;
                    default:
                        x[i] = 3; break;
                }
                ++count;
                if (count == 3) count = 0;
            }
            return x;
        }
        else
        {
            Console.WriteLine("非法变量数");
            return x;
        }
    }

    //计算收敛限
    static bool GetTol(int scope, double[] x)
    {
        List<double> f = new List<double>();
        double tol_cur = 0;
        double[] tol_den = GetFuncDen(scope, x);
        for (int i = 0; i < scope + 2; i++)
        {
            tol_cur += Math.Pow(tol_den[i], 2);
        }
        tol_cur = Math.Sqrt(tol_cur);
        fig.Add(tol_cur);
        return tol_cur < tol;
    }
    #endregion

    public static void Main(string[] args)
    {
        int[] count = new int[3];

        double[] x_ini1 = GetX(scope);
        double[] x11 = new double[scope + 2];
        for (int l = 0; l < iter; l++)
        {
            if (l == 0)
            {
                x11 = CalGrads(scope, x_ini1);
            }
            else
            {
                x11 = CalGrads(scope, x11);
            }
            if (GetTol(scope, x11))
            {
                count[0] = l; break;
            }
        }

        double[] x_ini2 = GetX(scope);
        double[] x22 = new double[scope + 2];
        for (int l = 0; l < iter; l++)
        {
            if (l == 0)
            {
                x22 = CalConjGrads(scope, x_ini2, l);
            }
            else
            {
                x22 = CalConjGrads(scope, x22, l);
            }
            if (GetTol(scope, x22))
            {
                count[1] = l; break;
            }
        }

        double[] x_ini3 = GetX(scope);
        double[] x33 = new double[scope + 2];
        for (int l = 0; l < iter; l++)
        {
            if (l == 0)
            {
                x33 = CalGradsBB(scope, x_ini3, l);
            }
            else
            {
                x33 = CalGradsBB(scope, x33, l);
            }
            if (GetTol(scope, x33))
            {
                count[2] = l; break;
            }
        }

        foreach (var item in fig)
        {
            Console.WriteLine(item);
        }
    }


    #region 梯度下降法
    //此函数用于迭代梯度下降法中的x
    static double[] CalGrads(int scope, double[] x)
    {
        double[] s = new double[scope + 2];
        s = GetFuncDen(scope, x);
        double[] xm = new double[scope + 2];
        double p = 0.2;
        double q = 0.5;
        for (int j = 0; j < 100; j++)
        {
            for (int i = 0; i < scope + 2; i++)
            {
                xm[i] = x[i] - q * s[i];
            }
            if (GetFuncNum(scope, x) - GetFuncNum(scope, xm) > 0.000000001)
            {
                break;
            }
            else
            {
                q = q * p;
            }
        }

        return xm;
    }
    #endregion

    #region 共轭梯度法
    static double g0 = 0;
    static double[] d0 = new double[scope + 2];

    //此函数用于迭代共轭梯度法中的x,选用的是Fletcher-Reeves公式
    static double[] CalConjGrads(int scope, double[] x, int l)
    {
        double[] xm = new double[scope + 2];
        double[] g = new double[scope + 2];
        double[] d = new double[scope + 2];
        double[][] G = GetG(scope);
        for (int i = 0; i < scope + 2; i++)
        {
            for (int j = 0; j < scope + 2; j++)
            {
                g[i] += G[i][j] * x[j];
                if (l == 0)
                {
                    d[i] -= G[i][j] * x[j];
                }
            }
        }
        double p = 0; double q = 0; double[] q1 = new double[scope + 2];
        if (l != 0)
        {
            double b = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                p += g[i] * g[i];
            }
            b = p / g0; g0 = p;
            for (int i = 0; i < scope + 2; i++)
            {
                d[i] = -g[i] + b * d0[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    q1[i] += d[j] * G[j][i];
                }
                q += d[i] * q1[i];
            }
            p /= q;
            for (int i = 0; i < scope + 2; i++)
            {
                d0[i] = d[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                xm[i] = x[i] + p * d[i];
            }
        }
        else
        {
            for (int i = 0; i < scope + 2; i++)
            {
                p += g[i] * g[i];
                for (int j = 0; j < scope + 2; j++)
                {
                    q1[i] += d[j] * G[j][i];
                }
                q += d[i] * q1[i];
            }
            g0 = p;
            p /= q; 
            for (int i = 0; i < scope + 2; i++)
            {
                d0[i] = d[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                xm[i] = x[i] + p * d[i];
            }
        }
        return xm;
    }

    static double[][] GetG(int scope)
    {
        double[][] G = new double[scope + 2][];
        //构造矩阵G
        for (int i = 0; i < scope + 2; i++)
        {
            G[i] = new double[scope + 2];
            for (int j = 0; j < scope + 2; j++)
            {
                if (i == 0)
                {
                    G[i][0] = 6; G[i][1] = -2; G[i][2] = -2;
                }
                else if (i == 1)
                {
                    G[i][0] = -2; G[i][1] = 12; G[i][2] = -4; G[i][2] = -2;
                }
                else if (i == scope)
                {
                    G[i][scope - 2] = -2; G[i][scope - 1] = -4; G[i][scope] = 12; G[i][scope + 1] = -2;
                }
                else if (i == scope + 1)
                {
                    G[i][scope - 1] = -2; G[i][scope] = -2; G[i][scope + 1] = 6;
                }
                else
                {
                    G[i][i] = 18; G[i][i - 1] = -4; G[i][i - 2] = -2; G[i][i + 1] = -4; G[i][i + 2] = -2;
                }
            }
        }
        return G;
    }
    #endregion

    #region 梯度BB法
    static double[] df1 = new double[scope + 2];
    static double[] x1 = new double[scope + 2];

    //此函数用于迭代梯度BB法中的x,选取的是第一类BB算子
    static double[] CalGradsBB(int scope, double[] x, int l)
    {
        double[] xm = new double[scope + 2];
        double[] df = new double[scope + 2];
        if (l == 0)
        {
            df = GetFuncDen(scope, x);
            xm = CalGrads(scope, x);
            for (int i = 0; i < scope + 2; i++)
            {
                df1[i] = df[i];
            }
            df = GetFuncDen(scope, xm);
            double p = 0; double q = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                p += (xm[i] - x[i]) * (df[i] - df1[i]);
                q += Math.Pow(df[i] - df1[i], 2);
            }
            p /= q;
            for (int i = 0; i < scope + 2; i++)
            {
                x1[i] = xm[i];
                xm[i] = xm[i] - p * df[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                df1[i] = df[i];
            }
        }
        else
        {
            df = GetFuncDen(scope, x);
            double p = 0; double q = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                p += (x[i] - x1[i]) * (df[i] - df1[i]);
                q += Math.Pow(df[i] - df1[i], 2);
            }
            p /= q;
            for (int i = 0; i < scope + 2; i++)
            {
                x1[i] = x[i];
                xm[i] = x[i] - p * df[i];
            }
            for (int i = 0; i < scope + 2; i++)
            {
                df1[i] = df[i];
            }
        }
        return xm;
    }
    #endregion
}