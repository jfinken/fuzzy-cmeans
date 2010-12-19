(ns fuzzy-cmeans.core)

;------------------------------------------------------------------------------
; ClusterPoint definition using datatypes and protocols 
; - to get the cluster-index for example: (:cluster-index my-point)
; - to update it: (assoc pts 0 (assoc (pts 0) :cluster-index 0.12345))
;   where pts is a seq of ClusterPoints
;------------------------------------------------------------------------------
(defprotocol cluster-point
  "protocol for cluster-points"
  (coords-at [point indx])
  (dimension [point]))

(defrecord ClusterPoint [coords cluster-index]
  cluster-point
  (coords-at [this indx] (coords indx))
  (dimension [this] (count coords)))

(defn make-cluster-point
  "Factory function returning a ClusterPoint object.  defrecord
  defines true Java classes, thus this function to prevent :import"
  [coords cluster-index]
  (ClusterPoint. coords cluster-index))

(defn update-cluster-index
  "Return a new ClusterPoint with the cluster-index assigned to val"
  [cpoint val]
  (assoc cpoint :cluster-index val))

(defn update-coord
  "Returns a new ClusterPoint with coordinate at k updated to val"
  [point k val]
  (assoc point :coords (assoc (:coords point) k val)))

; Helper fns for the U matrix
(defn val-colrow-of-vec
  [col row num-cols vec]
  (vec (+ col (* num-cols row))))
  
(defn alter-colrow-of-ref-vec
  [col row num-cols ref-vec value]
  (dosync 
    (ref-set ref-vec 
      (assoc @ref-vec (+ col (* num-cols row)) value))))

;------------------------------------------------------------------------------
; initialization routines 
;------------------------------------------------------------------------------
(def data-points (ref []))
(def clusters (ref []))
(def U (ref []))
(def fuzzy (ref -1))
(def eps (ref (Math/pow 10 -9))) ; default algorithm precision

; accessor routines, an iterface
(defn get-data-points
  "Returns the vector of ClusterPoints currently being operated on
  by fuzzy-cmeans"
  []
  @data-points)

(defn get-cluster-centroids
  "Returns the vector of cluster centroids currently being used in
  fuzzy-cmeans"
  []
  @clusters)

(defn U-at
  "Value of the U matrix at i j"
  [row col]
  (val-colrow-of-vec col row (count @clusters) @U))

; sum a vector
(defn sum-vec
  "Simply sum the contents of a vec"
  [vec]
  (apply + vec))

; euclidean distance in only 2-d
(defn euclid-dist-2d-or-eps
  "Return the 2D euclidean distance of two cluster points"
  [cp1 cp2]
  (let [diff
        (Math/sqrt (+ (Math/pow (- (coords-at cp1 0) (coords-at cp2 0)) 2)
                      (Math/pow (- (coords-at cp1 1) (coords-at cp2 1)) 2)))]
    (if (zero? diff)
      @eps
      diff)))
       
(defn init-U-matrix-1
  "Part 1 of 3 in the initialization of the U matrix"
  [cpoint i]
  (doseq [j (range (count @clusters))]
    (alter-colrow-of-ref-vec j i 
      (count @clusters) U (euclid-dist-2d-or-eps cpoint (@clusters j)))
    ))

(defn init-U-matrix-2
  "Part 2 of 3 in the initialization of the U matrix"
  [cpoint i]
 (let [curr-sum (sum-vec @U)]
   (doseq [j (range (count @clusters))]
     (alter-colrow-of-ref-vec j i
       (count @clusters) U 
       ; value to stuff into U
       (/ 1.0 (Math/pow (/ (U-at i j) curr-sum) (/ 2.0 (- @fuzzy 1.0)))))
     )))

(defn init-U-matrix-3
  "Part 3 of 3 in the initialization of the U matrix"
  [cpoint i]
  (let [curr-sum (sum-vec @U)]
    (doseq [j (range (count @clusters))]
      (alter-colrow-of-ref-vec j i
        (count @clusters) U
        ; value to stuff into U
        (/ (U-at i j) curr-sum)))
    ))

(defn recalc-cluster-index
  "With an initialized U matrix, update the cluster index of the given point"
  [cpoint i]
  (let [max (ref -1)]
    (doseq [j (range (count @clusters))]
      (if (< @max (U-at i j))
        ; do two things: set max to Uij and update cluster index
        (do
          ; alter max
          (dosync (ref-set max (U-at i j)))
          ; update cluster index for the point
          (if (== @max 0.5)
            (dosync (ref-set data-points (assoc @data-points i (update-cluster-index cpoint 0.5))))
            (dosync (ref-set data-points (assoc @data-points i (update-cluster-index cpoint j))))))
        nil))))

(defn calculate-cluster-indices
  []
  (doseq [i (range (count @data-points))]
    (recalc-cluster-index (@data-points i) i))
  )

(defn init-cmeans
  "Init the algorithm with a list of cluster-points, a list of cluster-points
  representing the initial number of clusters and their centroids, an initial
  fuzzy value, and an accuracy value (i.e. (Math/pow 10 -9))
  
  A large fuzzy value results in smaller memberships Uij and hence fuzzier clusters"
  [in-points in-clusters in-fuzzy, in-eps]
  ; stash
  (dosync (ref-set data-points in-points))
  (dosync (ref-set clusters in-clusters))
  (dosync (ref-set fuzzy in-fuzzy))
  (dosync (ref-set eps in-eps))
  
  ; loop over the points and clusters to create the initial U matrix
  (doseq [i (range (count @data-points))]
    (init-U-matrix-1 (@data-points i) i)
    (init-U-matrix-2 (@data-points i) i)
    (init-U-matrix-3 (@data-points i) i))
  
  ; recalculate cluster indices
  (calculate-cluster-indices)
  )

;------------------------------------------------------------------------------
; Calculate objective function routines
;------------------------------------------------------------------------------
(defn euler-distance
  "Return the distance between the point and cluster centroid"
  [point centroid]
  (loop [sum 0.0 i 0]
    (if (> i (- (dimension point) 1))
      (Math/sqrt sum)
      (recur 
        (+ sum (Math/pow (- ((:coords point) i) ((:coords centroid) i)) 2))
        (inc i)))))

(defn summed-cluster-distance
  "Helper routine for calculate-objective-function"
  [cpoint i]
  (loop [sum 0.0 j 0]
    (if (> j (- (count @clusters) 1))
      sum
      (recur
        (+ sum (* (Math/pow (U-at i j) @fuzzy)
                  (Math/pow (euler-distance cpoint (@clusters j)) 2)))
        (inc j)))))

(defn calculate-objective-function
  "Each iteration is based on minimizing this function.  The function represents
   the distance from any given data point to a cluster center weighted by that 
   data point's membership grade (U matrix)"
  []
  (loop [jk 0.0 i 0]
    (if (> i (- (count @data-points) 1))
      jk
      (recur
        (+ jk (summed-cluster-distance (@data-points i) i))
        (inc i)))))

;------------------------------------------------------------------------------
; Re-calculate cluster centers
;------------------------------------------------------------------------------
(defn calc-cluster-center-numer
  [centroid j]
  (let [uC (ref (vec (repeat (dimension centroid) 0)))]
    (doseq [i (range (count @data-points))]
      (let [uu (Math/pow (U-at i j) @fuzzy)]
        (doseq [k (range (dimension centroid))]
          (dosync (ref-set uC (assoc @uC k (+ (@uC k) (* uu ((:coords centroid) k)))))))))
    uC))

(defn calc-cluster-center-denom
  [centroid j]
  (loop [l 0 i 0]
    (if (> i (- (count @data-points) 1))
      l
      (recur
        (+ l (Math/pow (U-at i j) @fuzzy))
        (inc i)))))

(defn calculate-cluster-centers
  "Do just that, calculate new cluster centers"
  []
  (doseq [j (range (count @clusters))]
    (let [centroid (@clusters j)]
      (let [uC (calc-cluster-center-numer centroid j)
            l  (calc-cluster-center-denom centroid j)]
        (doseq [k (range (dimension centroid))]
          ; update the coordinates for the centroid
          (dosync (ref-set clusters (assoc @clusters j (update-coord centroid k
                                                        (/ (uC k) l)))))
          )))))
          
;------------------------------------------------------------------------------
; Perform one step of the algorithm 
;------------------------------------------------------------------------------
(defn update
  []
  (doseq [c (range (count @clusters))]
    (doseq [h (range (count @data-points))]
      (let [t1 (euler-distance (@data-points h) (@clusters c))]
        (let [top (if (> t1 1.0) t1 @eps) sumTerms (ref 0.0)]
          (doseq [ck (range (count @clusters))]
            (let [dist1 (euler-distance (@data-points h) (@clusters ck))]
              (let [dist (if (> dist1 1.0) dist1 @eps)]
                ; mutate sumTerms, as in sumTerms += blah
                (dosync (ref-set sumTerms (+ @sumTerms
                                             (Math/pow (/ top dist) (/ 2 (- @fuzzy 1)))))))))
          ; now update the U matrix
          (alter-colrow-of-ref-vec c h (count @clusters) U
            ; value to stuff into U
            (/ 1.0 @sumTerms))))))
  ; re-calculate the cluster centroids
  (calculate-cluster-centers))
  
;------------------------------------------------------------------------------
; Run!
;------------------------------------------------------------------------------
(defn run-iter
  "Perform n complete iterations of fuzzy-cmeans."
  [num-iters]
  (doseq [i (range num-iters)]
    (let [j (calculate-objective-function)]
      (calculate-cluster-centers)
      (update)
      (let [jnext (calculate-objective-function)]
        (println "iter" i "of fuzzy-cmeans! (j - jnext):" (Math/abs (- j jnext)))))))

(defn run
  "Perform a complete run of fuzzy-cmeans until the desired accuracy is
  achieved."
  []
  (loop [iter 0 accur 1000]
    (if (< accur @eps)
      accur
      (recur 
        (inc iter)
        (do 
          (println "iteration" iter "of fuzzy cmeans...")
          (let [j (calculate-objective-function)]
            (calculate-cluster-centers)
            (update)
            (let [jnext (calculate-objective-function)]
              (Math/abs (- j jnext)))))))))
;------------------------------------------------------------------------------
; Debug code
;------------------------------------------------------------------------------
(defn print-points
  []
  (doseq [i (range (count @data-points))]
    (println "point" i "cluster-index" (:cluster-index (@data-points i)) "location:" (:coords (@data-points i)))))

(defn print-clusters
  []
  (doseq [i (range (count @clusters))]
    (println "cluster" i "location: " (:coords (@clusters i)))))

(defn print-U-matrix
  []
  (doseq [i (range (count @data-points))]
    (doseq [j (range (count @clusters))]
      (println "i:"i "j:"j "Uij:" (U-at i j)))))

