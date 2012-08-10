(ns fuzzy-cmeans.test.core
    (:require [fuzzy-cmeans.core :as fuzzy] :reload)
    (:use [clojure.contrib.string :only (blank?)])
    (:use [incanter core stats charts])
    (:use [clojure.test]))

;(deftest replace-me ;; FIXME: write
;  (is false "No tests have been written."))

;------------------------------------------------------------------------------
; code to read a CSV of data points - currently unused
;------------------------------------------------------------------------------
(defn read-points-file
  [path]
  (with-open [rdr (clojure.java.io/reader path)]
    (reduce conj [] (line-seq rdr))))

(defn parse-point-str
  "Given a comma-separated point as a string, return
  a vector of doubles"
  [pt-str]
  (let [temp (.split pt-str ",")]
   (loop [side nil
          ret-vec (vector)
          more (seq temp)]
     (if (not (blank? (first more)))
       (recur 
         (println "first more:" (first more))
         (conj ret-vec (Double/parseDouble (first more)))
         (next more))
       ret-vec))))

(defn parse-points-file
  "Given a file where each line is a comma-sep list of dimensional
  values, construct a ClusterPoint for each"
  [path]
  (let [raw-pts-vec (read-points-file path)]
    (loop [ret-pts (vector)
           more (seq raw-pts-vec)]
      (if (boolean (seq (first more)))
        (recur (conj ret-pts (fuzzy/make-cluster-point (parse-point-str (first more)) -1))
               (next more))
        ret-pts))))

;------------------------------------------------------------------------------
; Client code
;------------------------------------------------------------------------------
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
                (fuzzy/make-cluster-point
                 (vector
                   ( + (mod (rand Integer/MAX_VALUE) (inc (- xmax xmin))) xmin) ;x
                   ( + (mod (rand Integer/MAX_VALUE) (inc (- ymax ymin))) ymin));y
                 -1))))))

(defn main
  []
  (let [in-fuzzy 2.0
        eps (Math/pow 10 -8) 
        xmin 1 xmax 500
        ymin 1 ymax 500
        num-clusters 5 
        num-points 3000]
    (let [data-pts (gen-cluster-points num-points xmin xmax ymin ymax)
          centroids (gen-cluster-points num-clusters xmin xmax ymin ymax)]
      ; init fuzzy-cmeans
      (fuzzy/init-cmeans data-pts centroids in-fuzzy eps)
      ; print clusters just for fun..
      (fuzzy/print-clusters)
      ; run 
      (fuzzy/run)
      )))
;------------------------------------------------------------------------------
; Visualize with Incanter 
;------------------------------------------------------------------------------
(defn third
  [more]
  (first (next (next more))))

(defn dump-points-for-incanter
  []
    (loop [i 0 ret (list)]
      (if (> i (- (count (fuzzy/get-data-points)) 1))
        ret
        (recur
          (inc i)
          (cons (conj (:coords ((fuzzy/get-data-points) i))
                      (:cluster-index ((fuzzy/get-data-points) i)))
                ret)))))

(defn dump-clusters-for-incanter
  []
    (loop [i 0 ret (list)]
      (if (> i (- (count (fuzzy/get-cluster-centroids)) 1))
        ret
        (recur
          (inc i)
          (cons (conj (:coords ((fuzzy/get-cluster-centroids) i))
                      (:cluster-index ((fuzzy/get-cluster-centroids) i)))
                ret)))))

(comment
(defn dump-point-cluster-indices
  []
  (loop [i 0 ret (vector)]
    (if (> i (- (count @data-points) 1))
      ret
        (recur
          (inc i)
          (conj ret (:cluster-index (@data-points i)))))))
)

(defn get-scatter-title
  [data]
  (str "fuzzy c-means (" (count data) " data points)"))

; grab the data and visualize with a scatter-plot
(defn incanter-vis
  []
  (let [data (dump-points-for-incanter)
        center-data (dump-clusters-for-incanter)]
    (let[ X (map first data)
          Y (map second data)
          cs (map third data)
          center_xs (map first center-data)
          center_ys (map second center-data)]
      (let [ plot (doto (scatter-plot X Y :group-by cs :title (get-scatter-title data))
                        (add-points center_xs center_ys))]
        (view plot)))))

; generate points and run fuzzy cmeans
(main)
; visualize
(incanter-vis)

