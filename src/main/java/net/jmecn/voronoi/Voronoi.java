package net.jmecn.voronoi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

public final class Voronoi {

    Random seeder;

    public List<Center> voronoiCenterList = new ArrayList<Center>();// voroni图所有的中心点
    public Map<Site, Corner> voronoiCornerList = new HashMap<Site, Corner>();// voroni图所有的拐角点
    public List<Edge> voronoiEdgeList = new ArrayList<Edge>();// vironoi图所有边

    public Voronoi() {
        seeder = new Random();
    }

    public void InitVoroni(int siteCount, int width, int height) {
        List<DelaunayTriangle> allTriangle = new ArrayList<DelaunayTriangle>();// delaunay三角形集合
        List<Site> sitesP = new ArrayList<Site>();
        int seed = seeder.nextInt();
        Random rand = new Random(seed);
        List<Edge> trianglesEdgeList = new ArrayList<Edge>();// Delaunay三角形网所有边

        List<Edge> voronoiRayEdgeList = new ArrayList<Edge>();// voroni图外围射线边
        // 初始设定点数为20
        // 初始设定画布大小是500*400
        // 超级三角形顶点坐标为（250,0），（0,400），（500,400）
        // 点集区域为（125,200），（125,400），（375,200），（375,400），随便设置，只要满足点落在三角形区域中
        for (int i = 0; i < siteCount; i++) {
            double x = (float) (rand.nextDouble() * width);
            double y = (float) (rand.nextDouble() * height);
            Site site = new Site(x, y);
            sitesP.add(site);
            Center c = new Center();
            c.point = site;
            voronoiCenterList.add(c);
        }

        // 按点集坐标X值排序
        sitesP.sort(new SiteSorterXY());
        voronoiCenterList.sort(new CenterSorterXY());
        for (int i = 0; i < voronoiCenterList.size(); i++) {
            voronoiCenterList.get(i).index = i;
        }

        int relaxNum = 5;
        for (int r = 0; r < relaxNum; r++) {
            if (r > 0) {
                allTriangle.clear();
                trianglesEdgeList.clear();
                voronoiEdgeList.clear();
                voronoiRayEdgeList.clear();
                voronoiCornerList.clear();

                for (int t = 0; t < voronoiCenterList.size(); t++) {
                    boolean isBorder = false;
                    double sumx = 0, sumy = 0;
                    double num = voronoiCenterList.get(t).corners.size();
                    for (Entry<Site, Corner> k : voronoiCenterList.get(t).corners.entrySet()) {
                        if (k.getValue().point.x <= 20 || k.getValue().point.x >= width - 20
                                || k.getValue().point.y <= 20 || k.getValue().point.y >= height - 20) {
                            isBorder = true;
                            break;
                        }
                        sumx += k.getValue().point.x;
                        sumy += k.getValue().point.y;
                    }
                    Center c = new Center();
                    if (!isBorder) {
                        c.point = new Site(sumx / num, sumy / num);
                        voronoiCenterList.set(t, c);
                        sitesP.set(t, new Site(sumx / num, sumy / num));
                    } else {
                        c.point = voronoiCenterList.get(t).point;
                        c.border = true;
                        voronoiCenterList.set(t, c);
                        sitesP.set(t, c.point);
                    }
                }

            }

            // 将超级三角形的三点添加到三角形网中
            Site A = new Site(0, 0);
            Site B = new Site(width, 0);
            Site C = new Site(0, height);
            Site D = new Site(width, height);
            Center CA = new Center();
            CA.point = A;
            CA.index = -1;
            Center CB = new Center();
            CB.point = B;
            CB.index = -1;
            Center CC = new Center();
            CC.point = C;
            CC.index = -1;
            Center CD = new Center();
            CD.point = D;
            CD.index = -1;
            sitesP.add(A);
            sitesP.add(B);
            sitesP.add(C);
            sitesP.add(D);

            voronoiCenterList.add(CA);
            voronoiCenterList.add(CB);
            voronoiCenterList.add(CC);
            voronoiCenterList.add(CD);

            DelaunayTriangle dt = new DelaunayTriangle(A, B, C, CA, CB, CC);
            DelaunayTriangle dt2 = new DelaunayTriangle(D, B, C, CD, CB, CC);
            allTriangle.add(dt);
            allTriangle.add(dt2);
            // 构造Delaunay三角形网
            setDelaunayTriangle(allTriangle, sitesP, voronoiCenterList);

            sitesP.sort(new SiteSorterXY());
            voronoiCenterList.sort(new CenterSorterXY());
            for (int i = 0; i < voronoiCenterList.size(); i++) {
                voronoiCenterList.get(i).index = i;
            }
            //
            // 不要移除，这样就不用画Delaunay三角形网外围边的射线
            // 移除超级三角形
            // voroObject.remmoveTrianglesByOnePoint(allTriangle, A);
            // voroObject.remmoveTrianglesByOnePoint(allTriangle, B);
            // voroObject.remmoveTrianglesByOnePoint(allTriangle, C);

            // 返回Delaunay三角形网所有边
            trianglesEdgeList = returnEdgesofTriangleList(allTriangle);

            // 填充neighbor
            for (int i = 0; i < allTriangle.size(); i++) {
                DelaunayTriangle t = allTriangle.get(i);
                t.center1.neighbors.put(t.center2.index, t.center2);
                t.center1.neighbors.put(t.center3.index, t.center2);
                t.center2.neighbors.put(t.center1.index, t.center1);
                t.center2.neighbors.put(t.center3.index, t.center3);
                t.center3.neighbors.put(t.center1.index, t.center1);
                t.center3.neighbors.put(t.center2.index, t.center2);
            }
            // 获取所有Voronoi边
            voronoiEdgeList = returnVoronoiEdgesFromDelaunayTriangles(allTriangle, voronoiRayEdgeList,
                    voronoiCenterList, voronoiCornerList);

            for (Entry<Site, Corner> k : voronoiCornerList.entrySet()) {
                Corner c = k.getValue();
                if (c.point.x <= 20 || c.point.x >= width - 20 || c.point.y <= 20 || c.point.y >= height - 20) {
                    c.border = true;
                }

            }
        }
    }

    // 根据Delaunay三角形网构造Voronoi图的边
    public List<Edge> returnVoronoiEdgesFromDelaunayTriangles(List<DelaunayTriangle> allTriangle,
            List<Edge> voronoiRayEdgeList, List<Center> voronoiCenterList, Map<Site, Corner> voronoiCornerList) {
        List<Edge> voronoiEdgeList = new ArrayList<Edge>();
        // List<Edge> voronoiRayEdgeList = new ArrayList<Edge>();
        for (int i = 0; i < allTriangle.size(); i++) {
            List<Edge> neighborEdgeList = new ArrayList<Edge>();// 三角形邻接边集合
            for (int j = 0; j < allTriangle.size(); j++)// 为了找出邻接三角形数为2的三角形，即最外边的三角形，循环只能从0开始
            {
                if (j != i)// 不与自身比较
                {
                    Edge neighborEdge = findCommonEdge(allTriangle.get(i), allTriangle.get(j));
                    if (neighborEdge != null) {
                        neighborEdgeList.add(neighborEdge);
                        // 构造Voronoi边
                        Corner c1, c2;
                        if ((c1 = voronoiCornerList.get(allTriangle.get(i).centerPoint)) == null) {
                            c1 = new Corner();
                            c1.point = allTriangle.get(i).centerPoint;
                            c1.index = voronoiCornerList.size();
                            voronoiCornerList.put(allTriangle.get(i).centerPoint, c1);
                        }
                        if ((c2 = voronoiCornerList.get(allTriangle.get(j).centerPoint)) == null) {
                            c2 = new Corner();
                            c2.point = allTriangle.get(j).centerPoint;
                            c2.index = voronoiCornerList.size();
                            voronoiCornerList.put(allTriangle.get(j).centerPoint, c2);
                        }
                        Edge voronoiEdge = new Edge(allTriangle.get(i).centerPoint, allTriangle.get(j).centerPoint, c1,
                                c2);
                        if (!voronoiEdgeList.contains(voronoiEdge)) {
                            voronoiEdgeList.add(voronoiEdge);

                            // c1的touches添加
                            c1.touches.put(neighborEdge.ca.index, neighborEdge.ca);
                            c1.touches.put(neighborEdge.cb.index, neighborEdge.cb);

                            // c2的touches添加
                            c2.touches.put(neighborEdge.ca.index, neighborEdge.ca);
                            c2.touches.put(neighborEdge.cb.index, neighborEdge.cb);

                            // c1的protrudes添加
                            c1.protrudes.put(voronoiEdge, voronoiEdge);

                            // c2的protrudes添加
                            c2.protrudes.put(voronoiEdge, voronoiEdge);

                            // c1的adjacent添加
                            c1.adjacent.put(c2.index, c2);

                            // c2的adjacent添加
                            c2.adjacent.put(c1.index, c1);

                            // A的corners添加
                            neighborEdge.ca.corners.put(c1.point, c1);
                            neighborEdge.ca.corners.put(c2.point, c2);

                            // B的corners添加
                            neighborEdge.cb.corners.put(c1.point, c1);
                            neighborEdge.cb.corners.put(c2.point, c2);

                            // A的borders添加
                            neighborEdge.ca.borders.put(voronoiEdge, voronoiEdge);

                            // B的borders添加
                            neighborEdge.cb.borders.put(voronoiEdge, voronoiEdge);

                        }
                    }
                }
            }
            if (neighborEdgeList.size() == 2)// 表示此三角形是外围三角形，Voronoi边需要射线
            {
                Site midpoint;
                Edge rayEdge = null;
                Corner c1, c2 = null;
                Center a = null, b = null;

                if ((c1 = voronoiCornerList.get(allTriangle.get(i).centerPoint)) == null) {
                    c1 = new Corner();
                    c1.point = allTriangle.get(i).centerPoint;
                    c1.index = voronoiCornerList.size();
                    voronoiCornerList.put(allTriangle.get(i).centerPoint, c1);
                }
                // 找出最外边并寻找中点构造Voronoi射线边
                if (isPointOnEdge(neighborEdgeList.get(0), allTriangle.get(i).site1)
                        && isPointOnEdge(neighborEdgeList.get(1), allTriangle.get(i).site1)) {
                    midpoint = findMidPoint(allTriangle.get(i).site2, allTriangle.get(i).site3);
                    rayEdge = produceRayEdge(allTriangle.get(i).centerPoint, midpoint);// 产生较长的射线，原理实现还是线段画出的线
                    voronoiRayEdgeList.add(rayEdge);
                    a = allTriangle.get(i).center2;
                    b = allTriangle.get(i).center3;

                }
                if (isPointOnEdge(neighborEdgeList.get(0), allTriangle.get(i).site2)
                        && isPointOnEdge(neighborEdgeList.get(1), allTriangle.get(i).site2)) {
                    midpoint = findMidPoint(allTriangle.get(i).site1, allTriangle.get(i).site3);
                    rayEdge = produceRayEdge(allTriangle.get(i).centerPoint, midpoint);
                    voronoiRayEdgeList.add(rayEdge);
                    a = allTriangle.get(i).center1;
                    b = allTriangle.get(i).center3;

                }
                if (isPointOnEdge(neighborEdgeList.get(0), allTriangle.get(i).site3)
                        && isPointOnEdge(neighborEdgeList.get(1), allTriangle.get(i).site3)) {
                    midpoint = findMidPoint(allTriangle.get(i).site1, allTriangle.get(i).site2);
                    rayEdge = produceRayEdge(allTriangle.get(i).centerPoint, midpoint);
                    voronoiRayEdgeList.add(rayEdge);
                    a = allTriangle.get(i).center1;
                    b = allTriangle.get(i).center2;

                }

                if (rayEdge != null && (c2 = voronoiCornerList.get(rayEdge.b)) == null) {
                    c2 = new Corner();
                    c2.point = rayEdge.b;
                    c2.index = voronoiCornerList.size();

                    voronoiCornerList.put(rayEdge.b, c2);
                }

                // c2的touches添加
                c2.touches.put(a.index, a);
                c2.touches.put(b.index, b);
                Edge borderEdge = new Edge(a.point, b.point, a, b);
                // c1的protrudes添加
                c1.protrudes.put(borderEdge, borderEdge);
                // c2的protrudes添加
                c2.protrudes.put(borderEdge, borderEdge);

                // c1的adjacent添加
                c1.adjacent.put(c2.index, c2);
                // c2的adjacent添加
                c2.adjacent.put(c1.index, c1);

                // A的corners添加
                a.corners.put(c2.point, c2);
                // B的corners添加
                b.corners.put(c2.point, c2);

                // A的borders添加
                a.borders.put(borderEdge, borderEdge);
                // B的borders添加
                b.borders.put(borderEdge, borderEdge);
            }
        }
        return voronoiEdgeList;
    }

    // 根据三角形链表返回三角形所有的边
    public List<Edge> returnEdgesofTriangleList(List<DelaunayTriangle> allTriangle) {
        List<Edge> commonEdges = new ArrayList<Edge>();
        for (int i = 0; i < allTriangle.size(); i++) {
            Edge edge1 = new Edge(allTriangle.get(i).site1, allTriangle.get(i).site2, allTriangle.get(i).center1,
                    allTriangle.get(i).center2);
            Edge edge2 = new Edge(allTriangle.get(i).site1, allTriangle.get(i).site3, allTriangle.get(i).center1,
                    allTriangle.get(i).center3);
            Edge edge3 = new Edge(allTriangle.get(i).site2, allTriangle.get(i).site3, allTriangle.get(i).center2,
                    allTriangle.get(i).center3);
            if (!commonEdges.contains(edge1))
                commonEdges.add(edge1);
            if (!commonEdges.contains(edge2))
                commonEdges.add(edge2);
            if (!commonEdges.contains(edge3))
                commonEdges.add(edge3);
        }
        return commonEdges;
    }

    // 根据点集构造Delaunay三角形网
    public void setDelaunayTriangle(List<DelaunayTriangle> allTriangle, List<Site> sites, List<Center> centers) {

        for (int i = 0; i < sites.size(); i++) {
            List<DelaunayTriangle> tmpTriList = new ArrayList<DelaunayTriangle>();
            // 拷贝所有三角形
            for (int j = 0; j < allTriangle.size(); j++) {
                tmpTriList.add(allTriangle.get(j));
            }

            // 受影响的三角形链表
            List<DelaunayTriangle> influenedTriangles = new ArrayList<DelaunayTriangle>();
            // 新形成的三角形链表
            List<DelaunayTriangle> newTriangles = new ArrayList<DelaunayTriangle>();
            // 受影响三角形的公共边
            List<Edge> commonEdges = new ArrayList<Edge>();

            for (int j = 0; j < tmpTriList.size(); j++) {
                double lengthToCenter;// 该点到圆心距离
                lengthToCenter = distance2Point(tmpTriList.get(j).centerPoint, sites.get(i));
                if (lengthToCenter < tmpTriList.get(j).radius) {
                    influenedTriangles.add(tmpTriList.get(j));// 添加到受影响的三角形链表
                    allTriangle.remove(tmpTriList.get(j));// 移除当前三角形
                }
            }

            // 从受影响的三角形链表中，形成新的三角形链表
            for (int k = 0; k < influenedTriangles.size(); k++) {
                addNewDelaunayTriangle(newTriangles, influenedTriangles.get(k), sites.get(i), centers.get(i));
            }

            // 查找受影响三角形的公共边
            if (influenedTriangles.size() > 1) {
                commonEdges = findCommonEdges(influenedTriangles);
            }
            // 将受影响三角形中的公共边所在的新形成的三角形排除
            if (commonEdges.size() > 0) {
                remmoveTrianglesByEdges(newTriangles, commonEdges);
            }
            // 对新形成的三角形进行局部优化
            LOP(newTriangles);
            // 将优化后的新形成的三角形添加到三角形链表中
            for (int k = 0; k < newTriangles.size(); k++) {
                allTriangle.add(newTriangles.get(k));
            }
        }
    }

    // 移除所有边边所在的三角形
    public void remmoveTrianglesByEdges(List<DelaunayTriangle> allTriangles, List<Edge> edges) {

        List<DelaunayTriangle> tmpTriList = new ArrayList<DelaunayTriangle>();
        // 拷贝所有三角形
        for (int i = 0; i < allTriangles.size(); i++) {
            tmpTriList.add(allTriangles.get(i));
        }

        for (int i = 0; i < tmpTriList.size(); i++) {
            for (int j = 0; j < edges.size(); j++) {
                if (isEdgeOnTriangle(tmpTriList.get(i), edges.get(j))) {
                    allTriangles.remove(tmpTriList.get(i));
                }
            }
        }
    }

    // 移除一条边所在的三角形
    public void remmoveTrianglesByOneEdge(List<DelaunayTriangle> allTriangles, Edge edge) {
        List<DelaunayTriangle> tmpTriList = new ArrayList<DelaunayTriangle>();
        // 拷贝所有三角形
        for (int i = 0; i < allTriangles.size(); i++) {
            tmpTriList.add(allTriangles.get(i));
        }
        for (int i = 0; i < tmpTriList.size(); i++) {
            if (isEdgeOnTriangle(tmpTriList.get(i), edge))
                allTriangles.remove(tmpTriList.get(i));
        }
    }

    // 移除点所在的三角形
    public void remmoveTrianglesByOnePoint(List<DelaunayTriangle> allTriangles, Site site) {
        List<DelaunayTriangle> tmpTriList = new ArrayList<DelaunayTriangle>();
        // 拷贝所有三角形
        for (int i = 0; i < allTriangles.size(); i++) {
            tmpTriList.add(allTriangles.get(i));
        }
        for (int i = 0; i < tmpTriList.size(); i++) {
            if (isPointOnTriangle(tmpTriList.get(i), site))
                allTriangles.remove(tmpTriList.get(i));
        }
    }

    // 判断边是否属于三角形
    public boolean isEdgeOnTriangle(DelaunayTriangle triangel, Edge edge) {
        int samePointNum = 0;
        if (siteIsEqual(edge.a, triangel.site1) || siteIsEqual(edge.a, triangel.site2)
                || siteIsEqual(edge.a, triangel.site3))
            samePointNum++;
        if (siteIsEqual(edge.b, triangel.site1) || siteIsEqual(edge.b, triangel.site2)
                || siteIsEqual(edge.b, triangel.site3))
            samePointNum++;
        if (samePointNum == 2)
            return true;
        return false;
    }

    // 判断点是否属于三角形
    public boolean isPointOnTriangle(DelaunayTriangle triangle, Site site) {
        if (siteIsEqual(site, triangle.site1) || siteIsEqual(site, triangle.site2) || siteIsEqual(site, triangle.site3))
            return true;
        return false;
    }

    // 判断点是否在边上
    public boolean isPointOnEdge(Edge edge, Site site) {
        if (siteIsEqual(site, edge.a) || siteIsEqual(site, edge.b))
            return true;
        return false;
    }

    // 将点与受影响的三角形三点连接，形成新的三个三角形添加到三角形链中
    public void addNewDelaunayTriangle(List<DelaunayTriangle> allTriangles, DelaunayTriangle influenedTri, Site point,
            Center centerPoint) {
        allTriangles.add(new DelaunayTriangle(influenedTri.site1, influenedTri.site2, point, influenedTri.center1,
                influenedTri.center2, centerPoint));
        allTriangles.add(new DelaunayTriangle(influenedTri.site1, influenedTri.site3, point, influenedTri.center1,
                influenedTri.center3, centerPoint));
        allTriangles.add(new DelaunayTriangle(influenedTri.site2, influenedTri.site3, point, influenedTri.center2,
                influenedTri.center3, centerPoint));
    }

    // 对新形成的三角形进行局部优化
    public List<DelaunayTriangle> LOP(List<DelaunayTriangle> newTriList) {
        List<DelaunayTriangle> resultTriList = new ArrayList<DelaunayTriangle>();
        // 拷贝新形成的三角
        for (int i = 0; i < newTriList.size(); i++) {
            resultTriList.add(newTriList.get(i));
        }

        for (int i = 0; i < newTriList.size(); i++) {
            for (int j = i + 1; j < newTriList.size(); j++) {
                Edge commonEdge;// 需要调整对角线的的三角形的公共边
                Site anotherPoint = new Site();// 新对角线的另一点
                Center anotherCenter = new Center();
                if (isInCircle(newTriList.get(j), newTriList.get(i).site1))// 三角形点在外接圆内
                {
                    // 找出两个三角形的公共边
                    commonEdge = findCommonEdge(newTriList.get(i), newTriList.get(j));
                    if (commonEdge != null) {
                        // 移除需要调整的三角形
                        resultTriList.remove(newTriList.get(i));
                        resultTriList.remove(newTriList.get(j));
                        // 找出对角线的另一点
                        if (siteIsEqual(newTriList.get(j).site1, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site1, commonEdge.b) == false) {
                            anotherPoint = newTriList.get(j).site1;
                            anotherCenter = newTriList.get(j).center1;
                        }
                        if (siteIsEqual(newTriList.get(j).site2, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site2, commonEdge.b) == false) {
                            anotherPoint = newTriList.get(j).site2;
                            anotherCenter = newTriList.get(j).center2;
                        }
                        if (siteIsEqual(newTriList.get(j).site3, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site3, commonEdge.b) == false) {
                            anotherPoint = newTriList.get(j).site3;
                            anotherCenter = newTriList.get(j).center3;
                        }
                        // 形成两个新的三角形
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site1, anotherPoint, commonEdge.a,
                                newTriList.get(i).center1, anotherCenter, commonEdge.ca));
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site1, anotherPoint, commonEdge.b,
                                newTriList.get(i).center1, anotherCenter, commonEdge.cb));
                    }
                }

                if (isInCircle(newTriList.get(j), newTriList.get(i).site2))// 三角形点在外接圆内
                {
                    // 找出两个三角形的公共边
                    commonEdge = findCommonEdge(newTriList.get(i), newTriList.get(j));
                    if (commonEdge != null) {
                        // 移除需要调整的三角形
                        resultTriList.remove(newTriList.get(i));
                        resultTriList.remove(newTriList.get(j));
                        // 找出对角线的另一点
                        if (siteIsEqual(newTriList.get(j).site1, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site1, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site1;
                        if (siteIsEqual(newTriList.get(j).site2, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site2, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site2;
                        if (siteIsEqual(newTriList.get(j).site3, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site3, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site3;
                        // 形成两个新的三角形
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site2, anotherPoint, commonEdge.a,
                                newTriList.get(i).center2, anotherCenter, commonEdge.ca));
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site2, anotherPoint, commonEdge.a,
                                newTriList.get(i).center2, anotherCenter, commonEdge.cb));
                    }
                }

                if (isInCircle(newTriList.get(j), newTriList.get(i).site3))// 三角形点在外接圆内
                {
                    // 找出两个三角形的公共边
                    commonEdge = findCommonEdge(newTriList.get(i), newTriList.get(j));
                    if (commonEdge != null) {
                        // 移除需要调整的三角形
                        resultTriList.remove(newTriList.get(i));
                        resultTriList.remove(newTriList.get(j));
                        // 找出对角线的另一点
                        if (siteIsEqual(newTriList.get(j).site1, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site1, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site1;
                        if (siteIsEqual(newTriList.get(j).site2, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site2, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site2;
                        if (siteIsEqual(newTriList.get(j).site3, commonEdge.a) == false
                                && siteIsEqual(newTriList.get(j).site3, commonEdge.b) == false)
                            anotherPoint = newTriList.get(j).site3;
                        // 形成两个新的三角形
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site3, anotherPoint, commonEdge.a,
                                newTriList.get(i).center3, anotherCenter, commonEdge.ca));
                        resultTriList.add(new DelaunayTriangle(newTriList.get(i).site3, anotherPoint, commonEdge.a,
                                newTriList.get(i).center3, anotherCenter, commonEdge.cb));
                    }
                }
            }
        }
        newTriList = resultTriList;
        return resultTriList;// 返回优化后的新形成的三角形
    }

    // 找出受影响的三角形的公共边
    public List<Edge> findCommonEdges(List<DelaunayTriangle> influenedTriangles) {
        List<Edge> coomonEdges = new ArrayList<Edge>();
        Edge tmpEdge;
        for (int i = 0; i < influenedTriangles.size(); i++) {
            for (int j = i + 1; j < influenedTriangles.size(); j++) {
                tmpEdge = findCommonEdge(influenedTriangles.get(i), influenedTriangles.get(j));
                if (tmpEdge != null) {
                    coomonEdges.add(tmpEdge);
                }
            }
        }
        return coomonEdges;
    }

    // 找出两个三角形的公共边
    public Edge findCommonEdge(DelaunayTriangle chgTri1, DelaunayTriangle chgTri2) {
        Edge edge;
        List<Site> commonSites = new ArrayList<Site>();
        List<Center> commonCenters = new ArrayList<Center>();
        if (siteIsEqual(chgTri1.site1, chgTri2.site1) || siteIsEqual(chgTri1.site1, chgTri2.site2)
                || siteIsEqual(chgTri1.site1, chgTri2.site3)) {
            commonSites.add(chgTri1.site1);
            commonCenters.add(chgTri1.center1);
        }
        if (siteIsEqual(chgTri1.site2, chgTri2.site1) || siteIsEqual(chgTri1.site2, chgTri2.site2)
                || siteIsEqual(chgTri1.site2, chgTri2.site3)) {
            commonSites.add(chgTri1.site2);
            commonCenters.add(chgTri1.center2);
        }
        if (siteIsEqual(chgTri1.site3, chgTri2.site1) || siteIsEqual(chgTri1.site3, chgTri2.site2)
                || siteIsEqual(chgTri1.site3, chgTri2.site3)) {
            commonSites.add(chgTri1.site3);
            commonCenters.add(chgTri1.center3);
        }
        if (commonSites.size() == 2) {
            edge = new Edge(commonSites.get(0), commonSites.get(1), commonCenters.get(0), commonCenters.get(1));
            return edge;
        }
        return null;
    }

    // 判断两点是否相同
    public boolean siteIsEqual(Site a, Site b) {
        if (a.x == b.x && a.y == b.y)
            return true;
        return false;
    }

    // 找出亮点的中点
    public Site findMidPoint(Site a, Site b) {
        Site midpoint = new Site((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
        return midpoint;

    }

    // 判断插入点是否在三角形边上
    public Site[] isOnEdges(DelaunayTriangle triangle, Site site) {
        Site[] edges = new Site[2];
        Site a = triangle.site1;
        Site b = triangle.site2;
        Site c = triangle.site3;

        if ((site.y - a.y) * (site.x - b.x) == (site.y - b.y) * (site.x - a.x))// 点在ab边上
        {
            edges[0] = a;
            edges[1] = b;
        }

        if ((site.y - a.y) * (site.x - c.x) == (site.y - c.y) * (site.x - a.x))// 点在ac边上
        {
            edges[0] = a;
            edges[1] = c;
        }

        if ((site.y - b.y) * (site.x - c.x) == (site.y - c.y) * (site.x - b.x))// 点在bc边上
        {
            edges[0] = b;
            edges[1] = c;
        }
        return edges;
    }

    // 判断点是否在三角形外接圆的内部
    public boolean isInCircle(DelaunayTriangle triangle, Site site) {
        double lengthToCenter;// 该点到圆心距离
        lengthToCenter = distance2Point(triangle.centerPoint, site);
        if (lengthToCenter < triangle.radius) {
            return true;
        }
        return false;
    }

    // 根据两点求以第一个点为起点的射线边
    public Edge produceRayEdge(Site start, Site direction) {
        Site end = new Site();
        Edge longEdge;
        end.x = 100 * (direction.x - start.x) + start.x;// 找出射线方向的较大的x终点
        end.y = (direction.y - start.y) * (end.x - start.x) / (direction.x - start.x) + start.y;
        longEdge = new Edge(start, end, new Corner(), new Corner());
        return longEdge;
    }

    // 求两点之间距离
    public double distance2Point(Site p, Site p2) {
        double value = Math
                .sqrt(Math.abs(p.x - p2.x) * Math.abs(p.x - p2.x) + Math.abs(p.y - p2.y) * Math.abs(p.y - p2.y));
        return value;
    }

    // 求三角形的外接圆心
    public double circle_center(Site center, Site sites0, Site sites1, Site sites2) {
        double x1, x2, x3, y1, y2, y3;
        double x = 0;
        double y = 0;

        x1 = sites0.x;
        x2 = sites1.x;
        x3 = sites2.x;
        y1 = sites0.y;
        y2 = sites1.y;
        y3 = sites2.y;
        x = ((y2 - y1) * (y3 * y3 - y1 * y1 + x3 * x3 - x1 * x1) - (y3 - y1) * (y2 * y2 - y1 * y1 + x2 * x2 - x1 * x1))
                / (2 * (x3 - x1) * (y2 - y1) - 2 * ((x2 - x1) * (y3 - y1)));
        y = ((x2 - x1) * (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1) - (x3 - x1) * (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1))
                / (2 * (y3 - y1) * (x2 - x1) - 2 * ((y2 - y1) * (x3 - x1)));

        center.x = x;
        center.y = y;
        double radius = Math.sqrt(
                Math.abs(sites0.x - x) * Math.abs(sites0.x - x) + Math.abs(sites0.y - y) * Math.abs(sites0.y - y));
        return radius;
    }

}