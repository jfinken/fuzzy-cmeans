(ns fuzzy-cmeans.core)

;------------------------------------------------------------------------------
; generating/working with a 2D array
;------------------------------------------------------------------------------
(comment
(into-array (map double-array [[1 2] [3 4]]))

(defn double-array-2d [coll]
  (let [w (count coll)
        h (apply max (map count coll))
        arr (make-array Double/TYPE w h)]
    (doseq [x (range w)
            y (range h)]
      (aset arr x y (double (get-in coll [x y]))))
    arr))


;; too complicated, what about a 1-d vector and sugar to
;; map 2d indexing into it.
;; i.e. for a cols x rows, m x n matrix:
;; (i, j) => col + (row * num_cols) => i + m*j
(def u (ref [2 4 6 8 10 12 14 16 18 20 22 24]))
)


;=> (ij-of-vec 2 2 4 3 @u)
;22
;=> (ij-of-vec 1 2 4 3 @u)
;20
; assignment
; (assoc [1 2 3] 3 45)

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
;(def eps (Math/pow 10 -5)) ; algorithm precision
(def eps 0.0001)

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
      eps
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
  representing the initial number of clusters and their centroids, and an initial
  fuzzy value.
  
  A large fuzzy value results in smaller memberships Uij and hence fuzzier clusters"
  [in-points in-clusters in-fuzzy]
  ; stash
  (dosync (ref-set data-points in-points))
  (dosync (ref-set clusters in-clusters))
  (dosync (ref-set fuzzy in-fuzzy))
  
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
        (let [top (if (> t1 1.0) t1 eps) sumTerms (ref 0.0)]
          (doseq [ck (range (count @clusters))]
            (let [dist1 (euler-distance (@data-points h) (@clusters ck))]
              (let [dist (if (> dist1 1.0) dist1 eps)]
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
(defn run
  "Perform a complete run of fuzzy-cmeans until the desired accuracy is
  achieved."
  []
  (doseq [i (range 5)]
    (let [j (calculate-objective-function)]
      (calculate-cluster-centers)
      (update)
      (let [jnext (calculate-objective-function)]
        (println "iter" i "of fuzzy-cmeans! (j - jnext):" (Math/abs (- j jnext)))))))


;------------------------------------------------------------------------------
; Client and debug code
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

(comment
; Generate random points 
(defn gen-cluster-points
  "Return a vector of n cluster points with random coordinate values
   within the given range"
  [n xmin xmax ymin ymax]
  (loop [m n ret (vector)]
    (if (zero? m)
      ret
      (recur 
        (dec m)
        (conj ret
                (ClusterPoint. 
                 (vector 
                   ( + (mod (rand Integer/MAX_VALUE) (inc (- xmax xmin))) xmin) ;x
                   ( + (mod (rand Integer/MAX_VALUE) (inc (- ymax ymin))) ymin));y
                 -1))))))

; sample, right off the bat
(def in-fuzzy 2.0)
(def xmin 1)
(def xmax 500)
(def ymin 1)
(def ymax 500)
(def num-clusters 5)
(def pts (gen-cluster-points 1000 xmin xmax ymin ymax))
;(def pts (get-test-data-points))
(def centroids (gen-cluster-points num-clusters xmin xmax ymin ymax))
;(def centroids (get-test-clusters))
; init!
(init-cmeans pts centroids in-fuzzy)
;(print-points)
(print-clusters)
; iterate
(run)
;------------------------------------------------------------------------------
; Incanter viz
;------------------------------------------------------------------------------
(defn third
  [more]
  (first (next (next more))))

(defn dump-points-for-incanter
  []
    (loop [i 0 ret (list)]
      (if (> i (- (count @data-points) 1))
        ret
        (recur
          (inc i)
          (cons (conj (:coords (@data-points i))
                      (:cluster-index (@data-points i)))
                ret)))))

(defn dump-clusters-for-incanter
  []
    (loop [i 0 ret (list)]
      (if (> i (- (count @clusters) 1))
        ret
        (recur
          (inc i)
          (cons (conj (:coords (@clusters i))
                      (:cluster-index (@clusters i)))
                ret)))))

(defn dump-point-cluster-indices
  []
  (loop [i 0 ret (vector)]
    (if (> i (- (count @data-points) 1))
      ret
        (recur
          (inc i)
          (conj ret (:cluster-index (@data-points i)))))))

(def data (dump-points-for-incanter))
(def center-data (dump-clusters-for-incanter))
(def xs (map first data))
(def ys (map second data))
(def cs (map third data))
;(def cs (dump-point-cluster-indices))
(def center_xs (map first center-data))
(def center_ys (map second center-data))

; scatter-plot the data points
;(def cx (vector ((:coords (@clusters 0)) 0)))
;(def cy (vector ((:coords (@clusters 0)) 1)))
(def plot (doto (scatter-plot xs ys :group-by cs)
            (add-points center_xs center_ys)))
(view plot)

)

(comment
(defn get-test-data-points
  []
  (conj (vector)
          (ClusterPoint. (vector 253, 482) -1)
          (ClusterPoint. (vector 123, 420) -1)
          (ClusterPoint. (vector 178, 99) -1)
          (ClusterPoint. (vector 346, 475) -1)
           (ClusterPoint. (vector 98, 495) -1)
           (ClusterPoint. (vector 483, 69) -1)
           (ClusterPoint. (vector 408, 395) -1)
           (ClusterPoint. (vector 32, 397) -1)
           (ClusterPoint. (vector 448, 327) -1)
           (ClusterPoint. (vector 295, 480) -1)
           (ClusterPoint. (vector 2, 169) -1)
           (ClusterPoint. (vector 414, 365) -1)
           (ClusterPoint. (vector 13, 251) -1)
           (ClusterPoint. (vector 77, 31) -1)
           (ClusterPoint. (vector 485, 343) -1)
           (ClusterPoint. (vector 469, 203) -1)
           (ClusterPoint. (vector 375, 212) -1)
           (ClusterPoint. (vector 304, 153) -1)
           (ClusterPoint. (vector 369, 57) -1)
           (ClusterPoint. (vector 58, 398) -1)
           (ClusterPoint. (vector 387, 194) -1)
           (ClusterPoint. (vector 157, 72) -1)
           (ClusterPoint. (vector 483, 393) -1)
           (ClusterPoint. (vector 63, 68) -1)
           (ClusterPoint. (vector 188, 377) -1)
           (ClusterPoint. (vector 91, 251) -1)
           (ClusterPoint. (vector 444, 312) -1)
           (ClusterPoint. (vector 457, 84) -1)
           (ClusterPoint. (vector 68, 257) -1)
           (ClusterPoint. (vector 408, 427) -1)
           (ClusterPoint. (vector 23, 315) -1)
           (ClusterPoint. (vector 490, 254) -1)
           (ClusterPoint. (vector 27, 280) -1)
           (ClusterPoint. (vector 193, 196) -1)
           (ClusterPoint. (vector 92, 378) -1)
           (ClusterPoint. (vector 29, 391) -1)
           (ClusterPoint. (vector 269, 397) -1)
           (ClusterPoint. (vector 93, 308) -1)
           (ClusterPoint. (vector 12, 342) -1)
           (ClusterPoint. (vector 381, 119) -1)
           (ClusterPoint. (vector 188, 9) -1)
           (ClusterPoint. (vector 342, 109) -1)
           (ClusterPoint. (vector 253, 218) -1)
           (ClusterPoint. (vector 258, 64) -1)
           (ClusterPoint. (vector 255, 220) -1)
           (ClusterPoint. (vector 85, 113) -1)
           (ClusterPoint. (vector 149, 131) -1)
           (ClusterPoint. (vector 376, 73) -1)
           (ClusterPoint. (vector 203, 403) -1)
           (ClusterPoint. (vector 45, 203) -1)))

(defn get-test-clusters
  []
  (conj (vector) 
        (ClusterPoint. (vector 253, 482) -1)
        (ClusterPoint. (vector 123, 420) -1)
        (ClusterPoint. (vector 178, 99) -1)))
  )
