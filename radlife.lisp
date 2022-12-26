;;; (load (compile-file "radlife.lisp"))

;;; ToDo:
;;;   Add effect range parameter

;;; This is a version of life that is based on a model of radioactive
;;; decay and devay chains. Rather than being a 1 or 0, as usual cells
;;; are numbers representing the atomic weight of the cell. Each element
;;; (i.e., type, i.e., number) has a probability of decaying along a set
;;; of specific paths, which can move them up or down the scale, and a
;;; normal time period for that decay. We don't actually model this very
;;; precisely for the moment, but have the intention of improving it.

(defparameter *wsize* 30)
(defparameter *ncells* (* *wsize* *wsize*))
(defparameter *world* (make-array (list *wsize* *wsize*)))
	
(defparameter *x* nil) ;; Output file

;;; For the moment we just make up the transition table randomly for
;;; testing purposes.

;;; Models and initializaion

(defvar *avg-delta-record* nil)

(defparameter *max-aw* 100) ;; how many elements you're going to get

(defparameter *aw->decay-model* (make-hash-table :test #'equal))
  
(defstruct dm aw pdown ndown pup nup)

;;; This is a confusing parameter. At, say, 0.05 it seems like it's
;;; telling you that 5% of your elements will be stable, but (a) it
;;; interacts with *max-p*, and (b) because there are two chances to
;;; be stable, i.e., to be below this threshold, once for pup and
;;; another for pdown, you are going to get twice the number you
;;; expect. So, to get ~5 stables in 100 (that is, with *max-aw* and
;;; *max-p* = 100) you would need to set this to 0.025.

(defparameter *stability-threshold* 0.20)

(defvar *stables* nil)

;;; Top of the pup or pdown scale (100=100%=1.0). You probably don't
;;; ever want this to be 100 because then there's going to be some
;;; element that always radiates, and if it's at the top or bottom of
;;; the scale, because we non-linearly threshold there, you'll never
;;; land -- those will just keep radiating and never walk down the
;;; ladder (because there's no where to go!) (??? Another way to do
;;; this would be to always force the highest and lowest elements to
;;; be stable.)

(defparameter *max-p* 80) 

(defparameter *show-stables-as-Xs* nil)

(defparameter *ladder-limit* 3) ;; Range of nup/ndown

;;; Just for Vlod:

(defmacro I-want-to-loop-your-baby! (&body baby)
  `(loop for x below *wsize*
	 do (loop for y below *wsize*
		  do ,@baby)))

(defun init-decay-models ()
  (setf *stables* nil)
  (format t "Initializing decay models:~%")
  (clrhash *aw->decay-model*)
  (format *x* "aw	pdown	ndown	pup	nup~%")
  (loop for aw below *max-aw*
	as dm = (make-dm :aw aw
			 :pdown (/ (random *max-p*) 100.0) :ndown (- (1+ (random *ladder-limit*)))
			 :pup (/ (random *max-p*) 100.0) :nup (1+ (random *ladder-limit*)))
	do
	;; A very low value for any of these should make the element
	;; stable, but it's unlikely to get a low value on both, so we
	;; hack it in any such case. (I think that this model is sort
	;; of wrong, but it's a first try.)
	(if (or (< (dm-pdown dm) *stability-threshold*)
		(< (dm-pup dm) *stability-threshold*))
	    (progn (push aw *stables*)
	      (setf (dm-pdown dm) 0.0
		    (dm-pup dm) 0.0)))
	(format *x* "~a	~a	~a	~a	~a~%"
		(dm-aw dm) (dm-pdown dm) (dm-ndown dm) (dm-pup dm) (dm-nup dm))
	(setf (gethash aw *aw->decay-model*) dm))
  (terpri *x*)
  )

(defun init-world ()
  (I-want-to-loop-your-baby! (setf (aref *world* x y) (random *max-aw*)))
  (print-*world*))

;;; Updating. Here's the good part: As usual we have a scan phase and
;;; then an update phase, but on the scan phase what we do is decide
;;; whether the cell decays on the current step. If it does, then the
;;; delta (up or down) is counted into the target/update array (which
;;; starts as all zeros), and a direction is chosen for the ejection
;;; of the emited particle, and the emitted value is added/subtracted
;;; from the current cell, and added to the cell in the selected
;;; direction. (For the moment we don't wrap around.) Then, in phase 2
;;; all we do is apply the deltas to the world.

(defvar *deltas* (make-array (list *wsize* *wsize*)))

(defun run-step ()
  ;; Init the delta array
  (init-delta-array)
  ;; Do the wild thing
  (radiate)
  ;; Finally, update the world
  (update-world)
  )
  
(defun init-delta-array ()
  (I-want-to-loop-your-baby! (setf (aref *deltas* x y) 0)))

(defun radiate ()
  (I-want-to-loop-your-baby! 
	(let*
	    ((dm (gethash (aref *world* x y) *aw->decay-model*))
	     ;; (NNN ndown is already negated at creation time)
	     (pdown (dm-pdown dm))
	     (pup (dm-pup dm))
	     (ndown (dm-ndown dm))
	     (nup (dm-nup dm))
	     (p (/ (random 100) 100.0)))
	  ;; We need to randomize the order we do this in
	  ;; otherwise it'll always head in one direction. UUU
	  (if (zerop (random 2))
	      (if (< p pdown)
		  (send-from x y ndown)
		  (if (< p pup)
		      (send-from x y nup)))
	      (if (< p pup)
		  (send-from x y nup)
		  (if (< p pdown)
		      (send-from x y ndown)))))))

(defun send-from (x y n)
  (incf (aref *deltas* x y) n)
  ;; Direction will be (+-1/0 . +-1/0). (NNN: This can sometimes
  ;; self-absorb, which will end up being a no-op (at least regrarding
  ;; its own raditation) because it will go up and down the same
  ;; amount. This is probably physically plausible -- virtual
  ;; radiation?)
  (let ((newx (+ x (pic-direction)))
	(newy (+ y (pic-direction))))
    ;; Edge protection (the edges are assumed to be absorbers)
    (when (and (>= newx 0) (>= newy 0)
	       (< newx *wsize*) (< newy *wsize*))
      (incf (aref *deltas* newx newy) n))))

(defun pic-direction ()
  (nth (random 3) '(-1 0 1)))

(defun update-world ()
  ;;(print '-----------)
  ;;(print *deltas*)
  (let ((delta-sum 0))
    (I-want-to-loop-your-baby! 
     (let ((new-val (max 0 (min (1- *max-aw*) (+ (aref *world* x y) (aref *deltas* x y))))))
       (incf delta-sum (aref *deltas* x y))
       (setf (aref *world* x y) new-val))
     )
    (let ((avg (float (/ delta-sum *ncells*))))
      (push avg *avg-delta-record*)
      (format *x* "~a	~a	~a	~a~%"
	      delta-sum avg (reduce #'+ *avg-delta-record*)
	      (/ (nstables) (float *ncells*))
	      )
      (push (copy-*world*) *worlds*)
      )))

;;; X is retained for stables
(defparameter *cset* "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWYZ1234567890)(*&^%$#@!=-+ []\{}|;:'\",./<>`~?")

(defparameter *lcset* (length *cset*))

(defparameter *worlds* nil)

(defun print-*world* ()
  (format t "~%-----------------------~%")
  (I-want-to-loop-your-baby!
   (let ((aw (aref *world* x y)))
     (format t " ~a" (cond ((and *show-stables-as-Xs* (member aw *stables*)) "X")
			   ((>= aw *lcset*) "_")
			   (t (aref *cset* aw)))))
   (format t "~%")))

(defun copy-*world* ()
  (let ((w (make-array (list *wsize* *wsize*))))
    (I-want-to-loop-your-baby! (setf (aref w x y) (aref *world* x y)))
    w))

;;; ===========================
;;; Gif conser

(defvar *ct* nil)

(defun gif-*worlds* (runid)
  (let* ((color-count 256)
         (color-table (skippy::make-color-table))
         (data-stream (skippy::make-data-stream :color-table color-table
						:loopingp t
						:height *wsize*
						:width *wsize*)))
    ;;; The color table is set up so that the stable elements have
    ;;; colors where the rgb are the same as their values, and the rest
    ;;; are not. 
    (dotimes (aw color-count)
      (skippy::add-color (decide-color-for-element aw) color-table))
    (loop for world in (reverse *worlds*)
	  as image-data = (make-array *ncells* :element-type '(unsigned-byte 8))
	  do (I-want-to-loop-your-baby! (setf (aref image-data (+ y (* x *wsize*))) (aref world x y)))
	  (skippy::make-image :height *wsize*
			      :width *wsize*
			      :data-stream data-stream ;; This will auto-add to the data-stream
			      :top-position 0
			      :left-position 0
			      :image-data image-data
			      :delay-time 5))
    (skippy::output-data-stream data-stream (format nil "results/~a_rl.gif" runid))))

(defun decide-color-for-element (aw)
  (let ((dm (gethash aw *aw->decay-model*)))
    (if (null dm) 
	(skippy::rgb-color (random 256) (random 256) (random 256))
	(if (zerop (dm-pup dm)) ;; Stable -> Grey
	    (skippy::rgb-color aw aw aw)
	    (skippy::rgb-color (random 256) (random 256) (random 256))))))

;;; Test

(defun run (&optional (nsteps *ncells*) &aux (steps 0) (runid (get-universal-time)))
  (with-open-file
      (*x* (format nil "results/~a_rl.xls" runid) :direction :output :if-exists :supersede)
    (setf *avg-delta-record* nil)
    (setf *worlds* nil)
    (init-decay-models)
    (init-world)
    (format *x* "delta_sum	instant_avg	sum_of_avgs	fract_sable~%")
    (block stable 
      (loop for i below nsteps
	    do (run-step) (incf steps)
	    (let ((nstables (nstables)))
	      ;;(when (zerop (random 1)) (print (list 'stables nstables)))
	      (when (= *ncells* nstables)
		(print "All atoms are stable!")
		(return-from stable nil)))))
    (print (print `((steps ,steps) ;; This is fun!!!
		    ,@(loop for var in '(*wsize* *ladder-limit* *show-stables-as-Xs* *max-p*
					 *max-aw* *ncells* *stability-threshold* *stables*)
			    collect (list var (eval var)))
		    (avg-avg-deltas ,(float (/ (reduce #'+ *avg-delta-record*) (length *avg-delta-record*))))
		    ))
	   *x*)
    (print-*world*)
    (print (list 'steps steps))
    )
  (gif-*worlds* runid)
  (regen-vhtml)
  )

(defun regen-vhtml ()
  (with-open-file
      (o "results/v.html" :direction :output :if-exists :supersede)
    (loop for file in (directory "results/*.gif")
	  as name = (pathname-name file)
	  do (format o "
<p><hr><p>
<h3>~a</h3>
<br>~s<br>
<image src=~a.gif width=400 height=400>~%
"
		     name (ignore-errors (get-info-block name)) name))))

(defun get-info-block (rlname) ;; UUU
  (with-open-file (i (format nil "results/~a.xls" rlname))
    (loop for val = (read i nil nil)
	  until (null val)
	  when (and (listp val) (eq 'steps (caar val)))
	  do (return val))))

(defun nstables (&aux (sum 0))
  (I-want-to-loop-your-baby! (if (member (aref *world* x y) *stables*) (incf sum)))
  sum)

(run 1000)
