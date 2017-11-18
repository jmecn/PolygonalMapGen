package net.jmecn.voronoi;

import java.util.List;

public class DelaunayTriangle {
    Voronoi voronoi = new Voronoi();
    public Site site1, site2, site3;// 三角形三点
    public Center center1, center2, center3;
    public Site centerPoint;// 外界圆圆心
    public double radius;// 外接圆半径
    public List<DelaunayTriangle> adjoinTriangle;// 邻接三角形

    public DelaunayTriangle(Site site1, Site site2, Site site3,Center center1, Center center2, Center center3)
    {
        centerPoint = new Site();
        this.site1 = site1;
        this.site2 = site2;
        this.site3 = site3;
        this.center1 = center1;

        this.center2 = center2;

        this.center3 = center3;

        //构造外接圆圆心以及半径
        radius = voronoi.circle_center(centerPoint, site1, site2, site3);
    }
}
