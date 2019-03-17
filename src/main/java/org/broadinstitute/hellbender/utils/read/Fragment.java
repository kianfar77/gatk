package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * All available evidence coming from a single biological fragment.  Either one read or a read pair.
 */
public class Fragment implements Locatable {

    private final SimpleInterval interval;

    private final List<GATKRead> reads;


    public Fragment(final GATKRead read) {
        reads = Collections.singletonList(read);
        interval = new SimpleInterval(read.getAssignedContig(), Math.min(read.getStart(), read.getEnd()), Math.max(read.getStart(), read.getEnd()));
    }

    public Fragment(final Pair<GATKRead, GATKRead> pair) {
        reads = Arrays.asList(pair.getLeft(), pair.getRight());
        final int start = Math.min(pair.getLeft().getStart(), pair.getRight().getStart());
        final int end = Math.max(pair.getLeft().getEnd(), pair.getRight().getEnd());
        interval = new SimpleInterval(pair.getLeft().getAssignedContig(), Math.min(start, end), Math.max(start, end));
    }

    public List<GATKRead> getReads() {
        return reads;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
