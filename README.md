# fuzzy-cmeans


## Overview
Fuzzy c-means (FCM) is a method of clustering which allows one piece of data 
to belong to two or more clusters, rather than belonging completely to one cluster. 
Points on the edge of a cluster may be in the cluster to a lesser degree than points 
in the center of cluster.

Here is an implementation in clojure.

## Usage

This implementation requires initialization with the data points to be sorted, 
the initial cluster centroids, and fuzzy and accuracy values.  The fuzzy parameter, a
minimum of 2.0, affects the degree to which the data points are members of the clusters.  A 
large fuzzy value technically results in smaller membership values and hence fuzzier 
clusters.

The fuzzy-cmeans algorithm supports clustering in n-dimensions (greater than zero!). Here
data points and cluster centroids are of type ClusterPoint implemented with defrecord.
However, a factory function exists to avoid having to use :import.

The below code from fuzzy-cmeans.test.core generates n random 2D data points via the
factory method make-cluster-point:

    (defn gen-cluster-points
      [n xmin xmax ymin ymax]
      (loop [m n ret (vector)]
        (if (zero? m)
          ret
          (recur
            (dec m)
              (conj ret
                (fuzzy/make-cluster-point
                  (vector
                    ( + (mod (rand Integer/MAX_VALUE) (inc (- xmax xmin))) xmin) ;x
                    ( + (mod (rand Integer/MAX_VALUE) (inc (- ymax ymin))) ymin));y
                  -1))))))

Below is an client implementation initializing and running fuzzy-cmeans: 

    (defn main
      []
      (let [in-fuzzy 2.0
            eps (Math/pow 10 -10)
            xmin 1 xmax 500
            ymin 1 ymax 500
            num-clusters 5
            num-points 1000]
        (let [data-pts (gen-cluster-points num-points xmin xmax ymin ymax)
              centroids (gen-cluster-points num-clusters xmin xmax ymin ymax)]
            ; init fuzzy-cmeans
            (fuzzy/init-cmeans data-pts centroids in-fuzzy eps)
            ; print clusters just for fun..
            (fuzzy/print-clusters)
            ; run 
            (fuzzy/run)
        )))

At this point I use Incanter to visualize how the clustered data points 
shake out.  Again please see fuzzy-cmeans.test.core.clj in test/

## License

Copyright (C) 2010 josh finken 

Distributed under the Eclipse Public License, the same as Clojure.
