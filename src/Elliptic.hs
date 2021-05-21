{-# LANGUAGE FlexibleInstances, OverlappingInstances #-}

module Elliptic where

-- |Representation of a group from abstract algebra
class Monoid a => Group a where
  (-.) :: a -> a
  (.-.) :: a -> a -> a
  x .-. y = x .+. (-.) y
  (.+.) :: a -> a -> a
  (.+.) = mappend

-- |Represent the part of the projective plane that is useful
-- to study for elliptic curves.
-- This means only one point at infinity for X=0 is added,
-- not for all directions (cos,sin)(theta)
data ProjectivePlane a = A a a
                       | O
  deriving (Show, Eq)

-- |Elliptic Curve defined using Weierstraß normal form
-- y^2 = x^3+ax^2+bx+c
data Elliptic a = Elliptic a a a
  deriving Eq

instance (Show a) => Show (Elliptic a) where
  show (Elliptic a b c) = "y^2 = x^3 + "
                        <> show a
                        <> "x^2 + "
                        <> show b
                        <> "x + "
                        <> show c


-- |Given an elliptic curve and a position on the x-Axis
-- return the y-Values, if they exist
yComp :: (Floating a, Ord a) => Elliptic a -> a -> Maybe (Either (a,a) a)
yComp (Elliptic a b c) x
  | v >  0 = Just $ Left (sqrt v, - sqrt v)
  | v == 0 = Just $ Right 0
  | otherwise = Nothing
  where
    v = a*x**2 + b*x + c

-- |Store an elliptic curve together with a point in the
-- projective Plane
type Solution a = (Elliptic a, ProjectivePlane a)

check :: (Floating a, Eq a) => Solution a -> Bool
check (Elliptic a b c, A x y) = y**2==x**3+a*x**2+b*x+c
check (Elliptic a b c, O) = True

(°) :: Floating a => Solution a -> Solution a -> Solution a
(e1@(Elliptic a b c), A x1 y1) ° (e2, A x2 y2) = (e1,A x3 y3)
  where
    x3 = lambda**2-a-x1-x2
    y3 = lambda*x3+nu
    lambda = (y2-y1)/(x2-x1)
    nu = y1-lambda*x1
(e1,A x y) ° (e2, O) = (e1, A x (-y))
v ° w = w ° v

-- |Halbgruppenstruktur auf elliptischen Kurven
-- Zugegeben, das ist etwas dumm, weil wir hier Tupel von
-- elliptischen Kurven mit Punkten auf einem pseudo-projektiven Raum
-- addieren, jedoch sind meine Haskell-Skills einfach noch nicht gut
-- genug für einfachere Lösungen
instance {-# OVERLAPS #-} (Eq a, Floating a) => Semigroup (Solution a) where
  p@(e@(Elliptic a b c),A x y) <> q
    | p/=q = (p ° q) ° (e,O)
    | p==q = (e, A (dup x y) (-sign x y *sqrt (dup x y)))
      where
        dup x y = (x**4-2*b*x**2-3*c*x+b**2-4*a*c)/(4*x**3+4*a*x**2+4*b*x+4*c)
        sign x y = signum y
  p <> q = q <> p

instance {-# OVERLAPS #-} (Eq a, Floating a) => Monoid (Solution a) where
  mempty = (Elliptic 0 0 0, O)

instance {-# OVERLAPS #-} (Eq a, Floating a) => Group (Solution a) where
  (-.) (e,A x y) = (e, A x (-y))
  (-.) (e,O) = (e, O)
