# Data section

Here are the downloaded data files. They come from [http://wsa.roe.ac.uk/].

Login credentials are still proprietary for most of these projects (6 July 2020).

Here is the base SQL query for one of the projects; essentially, the only thing that needs to be changed between any of them are the WSERV# references.

select SOURCEID, MEANMJDOBS, s.RA, s.DEC, JMHPNT, JMHPNTERR, HMKPNT, HMKPNTERR, JAPERMAG3, JAPERMAG3ERR, HAPERMAG3, HAPERMAG3ERR, KAPERMAG3, KAPERMAG3ERR, JPPERRBITS, HPPERRBITS, KPPERRBITS, MERGEDCLASS, PSTAR from                                         
wserv5SourceXSynopticSourceBestMatch as b, wserv5SynopticMergeLog as l, wserv5SynopticSource as s  where b.synFrameSetID=s.synFrameSetID and b.synSeqNum=s.synSeqNum and b.synFrameSetID=l.synFrameSetID and   s.RA > 0  order by SOURCEID, MEANMJDOBS

download as FITS, compress as GZIP.