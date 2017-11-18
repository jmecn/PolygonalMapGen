package net.jmecn.voronoi;

import java.util.Comparator;

public class CenterSorterXY implements Comparator<Center>
{
    public int compare(Center p1, Center p2)
    {
        if (p1.point.x > p2.point.x) return 1;
        if (p1.point.x < p2.point.x) return -1;
        return 0;
    }
}