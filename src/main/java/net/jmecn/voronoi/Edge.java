package net.jmecn.voronoi;
public class Edge
{

    public Site a, b;
    public Center ca, cb;
    public Corner cora, corb;
    public int river;
    public Edge(Site a, Site b,Center ca,Center cb)
    {
        this.a = a;
        this.b = b;
        this.ca = ca;
        this.cb = cb;

    }

    public Edge(Site a, Site b, Corner ca, Corner cb)
    {
        this.a = a;
        this.b = b;
        this.cora = ca;
        this.corb = cb;

    }
}