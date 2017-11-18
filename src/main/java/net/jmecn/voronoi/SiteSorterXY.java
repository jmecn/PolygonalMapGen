package net.jmecn.voronoi;

import java.util.Comparator;

public class SiteSorterXY implements Comparator<Site> {

    @Override
    public int compare(Site p1, Site p2)
    {
        if ( p1.x > p2.x ) return 1;
        if (p1.x < p2.x) return -1;
        return 0;
    }
}